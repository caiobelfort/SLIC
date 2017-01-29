#include "SLIC3D.h"

#include <itkGradientImageFilter.h>

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

	auto step = cbrt(num_voxels / m_number_of_super_pixels);

	int x_offset = step / 2;
	int y_offset = step / 2;
	int z_offset = step / 2;

	m_cluster_centers.clear();

	

	for(auto i = 0; i < width; ++i)
	{
		auto x = i * step;
		if (x > width - 1) break;
		for(auto j = 0 ; j < height; ++j)
		{
			auto y = j * step;
			if (y > height - 1) break;

			for(auto k = 0; k < depth; ++k)
			{
				auto z = k * step;
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

bool SLIC3D::IsOutOfBounds(float x, float y, float z)
{
	auto sz = m_input->GetLargestPossibleRegion().GetSize();
	return !(x >= 0 && x < sz[0] && y >= 0 && y < sz[1] && z >= 0 && z < sz[2]);
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

					//Check bound
					if (IsOutOfBounds(x, y, z)) continue;

					IndexType idx;
					idx[0] = x;
					idx[1] = y;
					idx[2] = z;

					auto gval = gimage->GetPixel(idx);
					IndexType cindex;
					cindex[0] = centroid.x;
					cindex[1] = centroid.y;
					cindex[2] = centroid.z;

					auto cgval = gimage->GetPixel(cindex);

					if(gval.GetNorm() < cgval.GetNorm())
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

void SLIC3D::DoSLICO()
{
	assert(m_input != nullptr);
	assert(m_output != nullptr);
	assert(m_level_distance_image != nullptr);
	assert(m_spatial_distance_image != nullptr);
	assert(m_distance_image != nullptr);
}

void SLIC3D::EnforceLabelConnectivity()
{
}

void SLIC3D::DrawContoursAroundSegmentedImage()
{
}
