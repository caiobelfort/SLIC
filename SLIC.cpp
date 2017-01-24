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
	m_distance = DistanceImageType::New();
	m_distance->SetRegions(m_input->GetLargestPossibleRegion());
	m_distance->Allocate();
	m_distance->FillBuffer(std::numeric_limits<float>::max());
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
	ImageType::IndexType index;


	//Cloning original project here
	int n(0); int r(0);
	for(auto j = 0; j < image_size[1]; ++j)
	{
		auto y = j * step + y_offset;
		if (y > image_size[1] - 1) break;

		for(auto i = 0; i < image_size[0]; ++i)
		{
			int x = i * step + (x_offset + r & 0x1); //Hex grid
			if (x > image_size[0] - 1) break;
			

			index[0] = x;
			index[1] = y;
			auto l = m_input->GetPixel(index);
			m_cluster_centers.push_back({ l, index });

			n++;
		}
		r++;
	}
	
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
				auto x = centroid.index[0] + i;
				auto y = centroid.index[1] + j;

				//Check bound
				if(x >= 0 && x < sz[0] && y >= 0 && y < sz[1])
				{
					GradientFilter::OutputImageType::IndexType idx, centroid_idx;

					idx[0] = x;
					idx[1] = y;

					auto gradient_value = gradient_image->GetPixel(idx);
					auto centroid_gradient_value = gradient_image->GetPixel(centroid.index);
					
					if (gradient_value.GetNorm() < centroid_gradient_value.GetNorm()) //Gradient is lower
					{
						centroid.index = idx;
						centroid.pixel_value = m_input->GetPixel(idx);
					}
				}
			}
		}
	}
}

void SLIC::DoSLICO()
{
	

}
