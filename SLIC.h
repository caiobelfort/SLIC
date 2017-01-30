#pragma once

#include <itkImage.h>

//First attempt unsigned char 2D
class SLIC
{
public:
	using ImageType = itk::Image<short, 2>;
	using LabeledImageType = itk::Image<int, 2>;
	using DistanceImageType = itk::Image<float, 2>;
	using IndexType = ImageType::IndexType;


	struct CentroidType
	{
		double l; 
		double x;
		double y;
	};

	void SetInput(ImageType::Pointer img){	m_input = img; }
	void SetNumberOfSuperPixels(int n) { m_number_of_super_pixels = n; }
	void SetPertubSeeds(bool flag) { m_pertub_seeds = flag;  }

	LabeledImageType::Pointer GetOutput() const { return m_output; }
	ImageType::Pointer GetCountourImage() const { return m_contour_image; }

    std::vector<ImageType::Pointer> GetImages() const;


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

    void GroupSimilarSuperpixels();
	
	//Computes the gradient of pixels and move centroids to lower gradient position in 3x3 neighborhood
	void PertubSeeds();


	void DoSLICO();

	void EnforceLabelConnectivity(std::vector<double> maxc, double inv);

	void DrawContoursAroundSegmentedImage();


};
