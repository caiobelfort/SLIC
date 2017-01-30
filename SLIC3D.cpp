#include "SLIC3D.h"

#include <itkGradientImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <iostream>
auto MakeIndex = [](auto x, auto y, auto z)
{
	SLIC3D::IndexType idx;
	idx[0] = x;
	idx[1] = y;
	idx[2] = z;

	return idx;
};

void SLIC3D::Update()
{

	InitCentroids();
	PertubSeeds();
	InitLabeledImage();
	InitDistanceImage();

	DoSLICO();
	EnforceLabelConnectivity();
	DrawContoursAroundSegmentedImage();
}

void SLIC3D::InitLabeledImage()
{
	m_output = LabeledImageType::New();
	m_output->SetRegions(m_input->GetLargestPossibleRegion());
	m_output->Allocate();
	m_output->FillBuffer(-1);
}

void SLIC3D::InitDistanceImage()
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

void SLIC3D::InitCentroids()
{
	auto num_voxels = m_input->GetLargestPossibleRegion().GetNumberOfPixels();
	auto image_size = m_input->GetLargestPossibleRegion().GetSize();
	auto width = image_size[0];
	auto height = image_size[1];
	auto depth = image_size[2];

	m_step = cbrt(num_voxels / m_number_of_super_pixels);

	int x_offset = m_step / 2;
	int y_offset = m_step / 2;
	int z_offset = m_step / 2;

	m_cluster_centers.clear();

	

	for(auto i = 0; i < width; ++i)
	{
		auto x = i * m_step;
		if (x > width - 1) break;
		for(auto j = 0 ; j < height; ++j)
		{
			auto y = j * m_step;
			if (y > height - 1) break;

			for(auto k = 0; k < depth; ++k)
			{
				auto z = k * m_step;
				if (z > depth - 1) break;


				IndexType index;
				index[0] = x;
				index[1] = y;
				index[2] = z;

				CentroidType centroid;
				centroid.l = m_input->GetPixel(index);
				centroid.x = x;
				centroid.y = y;
				centroid.z = z;
				m_cluster_centers.push_back(centroid);

			}
		}
	}
	m_number_of_super_pixels = m_cluster_centers.size();
}



void SLIC3D::PertubSeeds()
{
	using GradientFilter = itk::GradientImageFilter<ImageType>;
	auto filter = GradientFilter::New();

	filter->SetInput(m_input);
	filter->Update();

	auto gimage = filter->GetOutput();
	//For each initialized centroid, move it to lower gradient voxel in 3x3 neighborhood

	for(auto centroid: m_cluster_centers)
	{
		for(auto i = -1; i <= 1; ++i)
		{
			for(auto j = -1; j <= 1; ++j)
			{
				for(auto k = -1; k <= 1; ++k)
				{
					auto x = centroid.x + i;
					auto y = centroid.y + j;
					auto z = centroid.z + k;

					
			

					auto idx = MakeIndex(x, y, z);

					if (IsInsideImageRegion(idx))
					{
						auto gval = gimage->GetPixel(idx);
						auto cindex = MakeIndex(centroid.x, centroid.y, centroid.z);

						auto cgval = gimage->GetPixel(cindex);

						if (gval.GetNorm() < cgval.GetNorm())
						{
							centroid.x = x;
							centroid.y = y;
							centroid.z = z;
							centroid.l = m_input->GetPixel(idx);
						}
					}
				}
			}
		}
	}
}

void SLIC3D::DoSLICO()
{
	assert(m_input != nullptr);
	assert(m_output != nullptr);
	assert(m_level_distance_image != nullptr);
	assert(m_spatial_distance_image != nullptr);
	assert(m_distance_image != nullptr);

	auto image_size = m_input->GetLargestPossibleRegion().GetSize();
	auto image_region = m_input->GetLargestPossibleRegion();
	auto num_voxels = image_region.GetNumberOfPixels();

	m_step = cbrt(num_voxels / m_number_of_super_pixels);
	auto offset = m_step;
	if (m_step < 10) offset = m_step * 1.5;

	int width = image_size[0];
	int height = image_size[1];
	int depth = image_size[2];

	std::vector<int> clusters_size;
	std::vector<double> max_level(m_number_of_super_pixels, 10 * 10);
	std::vector<double> max_xyz(m_number_of_super_pixels, m_step * m_step);
	std::vector<CentroidType> centroid_sigma;

	auto invxyzwt = 1.0/(m_step * m_step);

	for(auto iterations = 0; iterations < 1; ++iterations)
	{
		std::cout << "Iteration: " << iterations << " | ";
		for(auto ck = 0; ck < m_cluster_centers.size(); ++ck)
		{	

			std::cout << "Superpixel: " << ck << std::endl;
			auto centroid = m_cluster_centers[ck];

			//Get search region
			

			auto xstart = std::max(0, static_cast<int>(centroid.x - offset));
			auto xend = std::min(width, static_cast<int>(centroid.x + offset));
			auto ystart = std::max(0, static_cast<int>(centroid.y - offset));
			auto yend = std::min(height, static_cast<int>(centroid.y + offset));
			auto zstart = std::max(0, static_cast<int>(centroid.z - offset));
			auto zend = std::min(depth, static_cast<int>(centroid.z + offset));

		

			//Iterate over search region
			for(auto x = xstart; x < xend; ++x)
			{
				for(auto y = ystart; y < yend; ++y)
				{
					for(auto z = zstart; z < zend; ++z)
					{
						auto index = MakeIndex(x, y, z);

						//Calculate distances
						m_level_distance_image->SetPixel(index, pow(m_input->GetPixel(index) - centroid.l, 2));
						m_spatial_distance_image->SetPixel(index, pow(x - centroid.x, 2) + pow(y - centroid.y, 2));
						auto distance = (m_level_distance_image->GetPixel(index) / max_level[ck]) +
							(m_spatial_distance_image->GetPixel(index) * invxyzwt);

						if (distance < m_distance_image->GetPixel(index))
						{
							m_distance_image->SetPixel(index, distance); //Change the distance
							m_output->SetPixel(index, ck); //Assign the label
						}

					}
				}
			}

			//Assign max color distance for a cluster
			if(iterations == 0)
			{
				max_level.assign(m_number_of_super_pixels, 1);
				max_xyz.assign(m_number_of_super_pixels, 1);
			}

			for(auto x = 0; x < width; ++x)
			{
				for(auto y = 0; y < height; ++y)
				{
					for(auto z = 0; z < depth; ++z)
					{
						auto index = MakeIndex(x, y, z);


						//pixel level distance
						if (max_level[m_output->GetPixel(index)] < m_level_distance_image->GetPixel(index))
							max_level[m_output->GetPixel(index)] = m_level_distance_image->GetPixel(index);

						//spatial distance
						if (max_xyz[m_output->GetPixel(index)] < m_distance_image->GetPixel(index))
							max_xyz[m_output->GetPixel(index)] = m_distance_image->GetPixel(index);
					}
				}
			}
			
			//Recalculate the centroids
			centroid_sigma.assign(m_number_of_super_pixels, { 0, 0, 0 });
			clusters_size.assign(m_number_of_super_pixels, 0);

			for (auto x = 0; x < width; ++x)
			{
				for (auto y = 0; y < height; ++y)
				{
					for (auto z = 0; z < depth; ++z)
					{
						auto index = MakeIndex(x, y, z);

						auto l = m_output->GetPixel(index);

						centroid_sigma[l].l += m_input->GetPixel(index);
						centroid_sigma[l].x += x;
						centroid_sigma[l].y += y;
						centroid_sigma[l].z += z;

						clusters_size[l]++;
					}
				}
			}

			for (auto k = 0; k < m_number_of_super_pixels; ++k)
			{
				if (clusters_size[k] <= 0) clusters_size[k] = 1;

				m_cluster_centers[k].l = centroid_sigma[k].l / clusters_size[k];
				m_cluster_centers[k].x = centroid_sigma[k].x / clusters_size[k];
				m_cluster_centers[k].y = centroid_sigma[k].y / clusters_size[k];
				m_cluster_centers[k].z = centroid_sigma[k].z / clusters_size[k];
			}
		}
	}


}

void SLIC3D::EnforceLabelConnectivity()
{
}

void SLIC3D::DrawContoursAroundSegmentedImage()
{	

	auto minMax = itk::MinimumMaximumImageCalculator<ImageType>::New();
	minMax->SetImage(m_input);
	minMax->Compute();

	auto max = minMax->GetMaximum();

	//Calculate gradients of labeled image
	auto gradient_filter = itk::GradientImageFilter<LabeledImageType>::New();
	gradient_filter->SetInput(m_output);
	gradient_filter->Update();

	//Generate a segmented image equal the input image
	m_contour_image = ImageType::New();
	m_contour_image->SetRegions(m_input->GetLargestPossibleRegion());
	m_contour_image->Allocate();
	m_contour_image->FillBuffer(0);

	auto image_size = m_input->GetLargestPossibleRegion().GetSize();

	auto gradient_image = gradient_filter->GetOutput();
	for (auto i = 0; i < image_size[0]; ++i)
	{
		for (auto j = 0; j < image_size[1]; ++j)
		{

			for (auto k = 0; k < image_size[2]; ++k)
			{
				auto index = MakeIndex(i, j, k);

				m_contour_image->SetPixel(index, m_input->GetPixel(index));

				if (gradient_image->GetPixel(index).GetNorm() > 0)
				{
					m_contour_image->SetPixel(index, max);
				}
			}
		}
	}
}
