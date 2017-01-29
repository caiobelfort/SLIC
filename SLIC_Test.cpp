
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "SLIC.h"
#include "SLIC3D.h"

int main(int argc, char **argv)
{
	if (argc < 2) return EXIT_FAILURE;

	using ImageType = itk::Image<short, 3>;

	auto reader = itk::ImageFileReader<ImageType>::New();
	reader->SetFileName(argv[1]);
	reader->Update();


	SLIC3D slic;
	slic.SetInput(reader->GetOutput());
	slic.SetNumberOfSuperPixels(200000);
	slic.Update();

	auto output = slic.GetCountourImage();

	return EXIT_SUCCESS;

}
