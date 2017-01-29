#pragma once
#include <itkImage.h>

class SLIC3D
{
public:
	using ImageType = itk::Image<short, 3>;
	using IndexType = ImageType::IndexType;
	using PixelType = ImageType::PixelType;
	using LabeledImageType = itk::Image<int, 3>;
	using DistanceImageType = itk::Image<float, 3>;

	struct CentroidType
	{
		float l;
		float x;
		float y;
		float z;
	};

	void SetInput(ImageType::Pointer img) { m_input = img; }
	void SetNumberOfSuperPixels(int n) { m_number_of_super_pixels = n; }
	void SetPertubSeeds(bool flag) { m_pertub_seeds = flag; }



	LabeledImageType::Pointer GetOutput() const { return m_output; }
	ImageType::Pointer GetCountourImage() const { return m_contour_image; }

	void Update();
	

private:

	ImageType::Pointer m_input = nullptr;
	ImageType::Pointer m_contour_image = nullptr;
	LabeledImageType::Pointer m_output = nullptr;
	DistanceImageType::Pointer m_level_distance_image = nullptr;
	DistanceImageType::Pointer m_spatial_distance_image = nullptr;
	DistanceImageType::Pointer m_distance_image = nullptr;

	int m_number_of_super_pixels = 200;
	bool m_pertub_seeds = false;
	std::vector<CentroidType> m_cluster_centers;


	void InitLabeledImage();
	void InitDistanceImage();
	void InitCentroids();

	bool IsOutOfBounds(float x, float y, float z);
	//Computes the gradient of pixels and move centroids to lower gradient position in 3x3 neighborhood
	void PertubSeeds();


	void DoSLICO();

	void EnforceLabelConnectivity();

	void DrawContoursAroundSegmentedImage();
};

