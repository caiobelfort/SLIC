#include "SLIC.h"

#include <itkGradientImageFilter.h>

void SLIC::Update()
{
	if (m_input == nullptr)
		throw std::runtime_error("Null input");

	InitCentroids();
	
	PertubSeeds();
	
	InitLabeledImage();
	
	InitDistanceImage();


	 
}

void SLIC::InitLabeledImage()
{
	m_output = LabeledImageType::New();
	m_output->SetRegions(m_input->GetLargestPossibleRegion());
	m_output->Allocate();
	m_output->FillBuffer(-1);
}

void SLIC::InitDistanceImage()
{
	m_distance_image = DistanceImageType::New();
	m_distance_image->SetRegions(m_input->GetLargestPossibleRegion());
	m_distance_image->Allocate();
	m_distance_image->FillBuffer(std::numeric_limits<float>::max());

	m_spatial_distance_image = DistanceImageType::New();
	m_spatial_distance_image->SetRegions(m_input->GetLargestPossibleRegion());
	m_spatial_distance_image->Allocate();
	m_spatial_distance_image->FillBuffer(std::numeric_limits<float>::max());

	m_level_distance_image = DistanceImageType::New();
	m_level_distance_image->SetRegions(m_input->GetLargestPossibleRegion());
	m_level_distance_image->Allocate();
	m_level_distance_image->FillBuffer(std::numeric_limits<float>::max());

}

void SLIC::InitCentroids()
{

	auto image_size = m_input->GetLargestPossibleRegion().GetSize();
	auto pixel_count = image_size[0] * image_size[1];

	//Calculate step
	auto step = sqrt(pixel_count / m_number_of_super_pixels);
	int x_offset = step / 2;
	int y_offset = step / 2;

	//Init cluster centers 
	m_cluster_centers.clear();

	//Cloning original project here
	int n(0); int r(0);
	for(auto j = 0; j < image_size[1]; ++j)
	{
		auto y = j * step + y_offset;
		if (y > image_size[1] - 1) break;

		for(auto i = 0; i < image_size[0]; ++i)
		{
			auto x = i * step + (x_offset + r & 0x1); //Hex grid
			if (x > image_size[0] - 1) break;
			
			IndexType index;
			index[0] = x;
			index[1] = y;

			CentroidType centroid;

			
			centroid.l = m_input->GetPixel(index);
			centroid.x = x;
			centroid.y = y;
			
			m_cluster_centers.push_back(centroid);

			n++;
		}
		r++;
	}
	
	m_number_of_super_pixels = m_cluster_centers.size();
}

void SLIC::PertubSeeds()
{	

	//Computes the gradient image
	using GradientFilter = itk::GradientImageFilter<ImageType>;
	auto gradient_image_filter = GradientFilter::New();
	gradient_image_filter->SetInput(m_input);
	gradient_image_filter->Update();
	auto gradient_image = gradient_image_filter->GetOutput();
	auto sz = m_input->GetLargestPossibleRegion().GetSize();

	//Now for each initialized centroid, search for lower gradient in 3x3 neighborhood and move to it
	for(auto centroid : m_cluster_centers)
	{
		for(auto j = -1 ; j <= 1; ++j)
		{	
			for(auto i = -1 ; i <= 1; ++i)
			{
				auto x = centroid.x + i;
				auto y = centroid.y + j;

				//Check bound
				if(x >= 0 && x < sz[0] && y >= 0 && y < sz[1])
				{
					GradientFilter::OutputImageType::IndexType idx;
					idx[0] = x;
					idx[1] = y;

					auto gradient_value = gradient_image->GetPixel(idx);
					IndexType index;
					index[0] = static_cast<int>(centroid.x);
					index[1] = static_cast<int>(centroid.y);
					auto centroid_gradient_value = gradient_image->GetPixel(index);
					
					if (gradient_value.GetNorm() < centroid_gradient_value.GetNorm()) //Gradient is lower
					{
						centroid.x = idx[0];
						centroid.y = idx[1];
						centroid.l = m_input->GetPixel(index);
					
					}
				}
			}
		}
	}
}

void SLIC::DoSLICO()
{
	assert(m_input != nullptr);
	assert(m_output != nullptr);
	assert(m_level_distance_image != nullptr);
	assert(m_spatial_distance_image != nullptr);
	assert(m_distance_image != nullptr);

	auto image_size = m_input->GetLargestPossibleRegion().GetSize();
	auto num_pixels = image_size[0] * image_size[1];
	auto step = sqrt(num_pixels / m_number_of_super_pixels);
	auto offset = step;
	if (step < 10) offset = step * 1.5;
	int width = image_size[0];
	int height = image_size[1];
	ImageType::IndexType idx, idy;


	std::vector<int> cluster_size;
	std::vector<double> max_level(m_number_of_super_pixels, 10 * 10);
	std::vector<double> max_xy(m_number_of_super_pixels, step*step);
	std::vector<CentroidType> centroid_sigma;

	auto invxywt = 1.0 / (step*step);

	for(auto n = 0; n < 10; ++n)
	{
		for(auto i = 0; i < m_cluster_centers.size(); ++i)
		{
			auto centroid = m_cluster_centers[i];
			//Get search region;
			idx[0] = std::max(0, static_cast<int>(centroid.x - offset));
			idx[1] = std::min(width, static_cast<int>(centroid.x + offset));
			idy[0] = std::max(0, static_cast<int>(centroid.y - offset));
			idy[1] = std::min(height, static_cast<int>(centroid.y + offset));


			//Iterate over the region
			for(auto x = idx[0]; x < idx[1]; ++x)
			{
				for(auto y = idy[0]; y < idy[1]; ++y)
				{
					ImageType::IndexType index;
					index[0] = x;
					index[1] = y;

				
					//Get the distances
					m_level_distance_image->SetPixel(index, pow(m_input->GetPixel(index) - centroid.l, 2));
					m_spatial_distance_image->SetPixel(index, pow(x - centroid.x, 2) + pow(y - centroid.y, 2));
					auto distance = (m_level_distance_image->GetPixel(index) / max_level[i]) + (m_spatial_distance_image->GetPixel(index) * invxywt);

					if(distance < m_distance_image->GetPixel(index))
					{
						m_distance_image->SetPixel(index, distance); //Change the distance
						m_output->SetPixel(index, i); //Assign the label
					}

				}
			}


		}

		// Assign the max color distance for a cluster
		if(n == 0)
		{
			max_level.assign(m_number_of_super_pixels, 1);
			max_xy.assign(m_number_of_super_pixels, 1);
		}
		for (auto x = 0; x < width; ++x)
		{
			for(auto y = 0; y < height; ++y )
			{	
				ImageType::IndexType index;
				index[0] = x;
				index[1] = y;
				
				//pixel level distance
				if (max_level[m_output->GetPixel(index)] < m_level_distance_image->GetPixel(index))
					max_level[m_output->GetPixel(index)] = m_level_distance_image->GetPixel(index);
				
				//spatial distance
				if (max_xy[m_output->GetPixel(index)] < m_distance_image->GetPixel(index))
					max_xy[m_output->GetPixel(index)] = m_distance_image->GetPixel(index);
			}
		}

		/**
		 * TODO
		 */
		//Centroids com sigma nao tao batendo!!!!!!!!

	}

	//repeat



	//for each centroid
	//for each pixel in 2S region of the centroid
	//
}
