
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "SLIC.h"

int main(int argc, char **argv)
{
	if (argc < 2) return EXIT_FAILURE;

	using ImageType = itk::Image<unsigned char, 2>;

	auto reader = itk::ImageFileReader<ImageType>::New();
	reader->SetFileName(argv[1]);
	reader->Update();


	SLIC slic;
	slic.SetInput(reader->GetOutput());
	slic.SetNumberOfSuperPixels(2000);
	slic.Update();

	auto output = slic.GetCountourImage();

	auto writer = itk::ImageFileWriter<SLIC::ImageType>::New();

	writer->SetFileName("segments.mha");
	writer->SetInput(output);
	writer->Update();

	return EXIT_SUCCESS;

}