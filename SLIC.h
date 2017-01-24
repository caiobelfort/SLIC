#pragma once

#include <itkImage.h>

//First attempt unsigned char 2D
class SLIC
{
public:
	using ImageType = itk::Image<unsigned char, 2>;
	using LabeledImageType = itk::Image<int, 2>;
	using DistanceImageType = itk::Image<float, 2>;
	

	struct Centroid
	{
		ImageType::PixelType pixel_value;
		ImageType::IndexType index;
	};

	void SetInput(ImageType::Pointer img){	m_input = img; }
	void SetNumberOfSuperPixels(int n) { m_number_of_super_pixels = n; }
	void SetPertubSeeds(bool flag) { m_pertub_seeds = flag;  }

	void Update();

private:
	ImageType::Pointer m_input = nullptr;
	LabeledImageType::Pointer m_output = nullptr;
	DistanceImageType::Pointer m_distance = nullptr;
	
	int m_number_of_super_pixels = 200;
	bool m_pertub_seeds = false;
	std::vector<Centroid> m_cluster_centers;

	void InitLabeledImage();
	void InitDistanceImage();
	void InitCentroids();
	
	//Computes the gradient of pixels and move centroids to lower gradient position in 3x3 neighborhood
	void PertubSeeds();


	void DoSLICO();

};