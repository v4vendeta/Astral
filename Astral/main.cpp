#include <iostream>
#include <string>
#include <glm/glm.hpp>
#include <boost/gil.hpp>
#include <boost/gil/extension/io/bmp.hpp>
#include <stdexcept>
#include <cstdint>
#include <filesystem>

int main(int argc, char **argv) {

	const uint32_t width{ 1280 }, height{ 720 };
	std::filesystem::path outpath{ "../output/ouput.bmp" };
	try {
		
		boost::gil::rgb8_image_t out{width,height};
		boost::gil::rgb8_pixel_t red(255, 255, 255);
		fill_pixels(boost::gil::view(out), red);
	
		boost::gil::write_view(outpath.c_str(), boost::gil::const_view(out), boost::gil::bmp_tag{});
		
	}
	catch (std::exception& e) {
		std::cout << e.what() << "\n";
	}
	return 0;

}