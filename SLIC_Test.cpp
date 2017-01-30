
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "SLIC.h"

using ImageType = itk::Image<short, 3>;
using SliceImageType = itk::Image<short, 2>;

SliceImageType::Pointer GetSlice(ImageType::Pointer img, int i)
{
    auto sz = img->GetRequestedRegion().GetSize();
    auto origin = img->GetOrigin();
    auto spacing = img->GetSpacing();
    SliceImageType::SizeType slice_sz;
    slice_sz[0] = sz[0];
    slice_sz[1] = sz[1];

    SliceImageType::IndexType slice_index;
    slice_index.Fill(0);

    SliceImageType::RegionType slice_region;
    slice_region.SetIndex(slice_index);
    slice_region.SetSize(slice_sz);

    auto slice = SliceImageType::New();

    slice->SetRegions(slice_region);
    slice->Allocate();

    for(auto x = 0 ; x < sz[0]; ++x)
    {
        for(auto y = 0; y < sz[1]; ++y)
        {
            ImageType::IndexType image_index;
            image_index[0] = x;
            image_index[1] = y;
            image_index[2] = i;

            SliceImageType::IndexType slice_index;

            slice_index[0] = x;
            slice_index[1] = y;

            slice->SetPixel(slice_index, img->GetPixel(image_index));

        }
    }

    SliceImageType::SpacingType slice_spacing;
    slice_spacing[0] = spacing[0];
    slice_spacing[1] = spacing[1];

    SliceImageType::PointType slice_origin;
    slice_origin[0] = origin[0];
    slice_origin[1] = origin[1];

    slice->SetOrigin(slice_origin);
    slice->SetSpacing(slice_spacing);
    return slice;

}

int main(int argc, char **argv)
{

	auto reader = itk::ImageFileReader<ImageType>::New();
	reader->SetFileName(argv[1]);
	reader->Update();


    auto img = reader->GetOutput();

    auto slice = GetSlice(img, 130);

	SLIC slic;
	slic.SetInput(slice);
	slic.SetNumberOfSuperPixels(1001);
	slic.Update();

	auto images = slic.GetImages();
    auto writer = itk::ImageFileWriter<SliceImageType>::New();
    for(auto i = 0 ; i < images.size(); ++i)
    {

	    writer->SetInput(images[i]);
        std::string name = std::to_string(i) + ".mha";
	    writer->SetFileName(name.c_str());
	    writer->Update();

    }
	return EXIT_SUCCESS;

}
