#include "Level.h"
#include <math.h>
#include <cstdint>
#include <cstdlib>
#include <complex>
#include <thread>
#include <mutex>
#include <amp.h>
#include <amp_math.h>
#include <array>
#include <vector>
#include <random>

#include "Complex_for_amp.h"

#define MAX_ITERATIONS 2000

using namespace concurrency;
using std::complex;


Level::Level(sf::RenderWindow* hwnd, Mouse* mus)
{
	window = hwnd;
	mouse = mus;
	
	// instead of running entire queery_amp i just ran first line of it
	//std::cout << accelerator::get_all().size() << std::endl;

	//compute_mandelbrot(-0.751085, -0.734975, 0.118378, 0.134488);
	//compute_mandelbrot_multicore(-2.0, 1.0, 1.125, -1.125);

	double* perlin = PerlinNoise::GeneratePerlin(); //output is values between 1 and -1
	//double* perlin = PerlinNoise::GeneratePerlinSerial();
	//double* perlin = PerlinNoise::GeneratePerlin_cpu_fromWiki();// used for reference

	perlin_image.create(PerlinNoise::WIDTH, PerlinNoise::HEIGHT, sf::Color::Magenta);

	for (int i = 0; i < PerlinNoise::HEIGHT; i++)
	{
		//std::cout << "\n new Row ";
		for (int j = 0; j < PerlinNoise::WIDTH; j++)
		{
			//	std::cout << perlin[i * PerlinNoise::WIDTH + j]+1 << " ";

			uint8_t R = 127 * (perlin[(i * PerlinNoise::WIDTH) + j] + 1);
			uint8_t G = 127 * (perlin[(i * PerlinNoise::WIDTH) + j] + 1);
			uint8_t B = 127 * (perlin[(i * PerlinNoise::WIDTH) + j] + 1);
			uint32_t color = (R << 24 | G << 16 | B << 8 | 255);

			perlin_image.setPixel(i, j, sf::Color(color));

		}

	}

	
	//perlin_tex.loadFromImage(blur(&perlin_image,13));
	perlin_tex.loadFromImage(perlin_image);


	//mandelbrot = new sf::Image;
	//mandelbrot->create(1024, 1024, sf::Color::Magenta);
	//compute_mandelbrot_gpu(-2.0, 1.0, 1.125, -1.125);

	//image palpatine
	//palpatine.loadFromFile("palpatine.jpg");
	//std::cout << palpatine.getSize().x << "x" << palpatine.getSize().y;
	//palpatine_tex.loadFromImage(blur(&palpatine,7));

	//rectangle
	//rect.setSize((sf::Vector2f)window->getSize());
	//rect.setSize(sf::Vector2f(PerlinNoise::WIDTH*2.5, PerlinNoise::WIDTH*2.5));
	rect.setTexture(&perlin_tex);
	//rect.setTexture(&mandelbrot_tex);
	//rect.setTexture(&palpatine_tex);
}

Level::~Level()
{
}

void Level::recalculate()
{
	//mouse position within a window
	sf::Vector2<double> mouse = (sf::Vector2<double>)(sf::Mouse::getPosition(*window));

	//align origins of both planes
	mouse.x -= (window->getSize().x /3.*2.);
	mouse.y -= (window->getSize().y/2.);

	//translate from one cordinate system to other

	mouse.x = (mouse.x / (double)window->getSize().x) * (3.);
	mouse.y = (mouse.y / (double)window->getSize().y) * (2.25);

	//mouse.x /= zoom;
	//mouse.y /= zoom;
	
	//debug
	//std::cout << mouse.x << " ";
	//std::cout << mouse.y << std::endl;

	//start coords
	left =    -2.0   /zoom + mouse.x;
	right =    1.0   /zoom + mouse.x;
	top =     -1.125 /zoom + mouse.y;
	bottom =   1.125 /zoom + mouse.y;

	//std::cout <<left<< " " << std::endl;

	compute_mandelbrot_gpu(left, right, top, bottom);
}

// handle user input
void Level::handleInput()
{
	/*if (sf::Mouse::isButtonPressed(sf::Mouse::Button::Left)) {
		if (!mouseDown)
		{
			mouseDown = true;
			offset = (sf::Vector2<double>)sf::Mouse::getPosition(*window);
		}
		

		if (mouseDown && !sf::Mouse::isButtonPressed(sf::Mouse::Button::Left)) {
			mouseDown = false;
		}

		recalculate();
	}*/

	//if (sf::Keyboard::isKeyPressed(sf::Keyboard::Z)) { zoom *= 1.2; recalculate();}
	//if (sf::Keyboard::isKeyPressed(sf::Keyboard::X) && zoom > 1) {	zoom /= 1.2; recalculate();	}
}

// Update game objects
void Level::update()
{
	//compute_mandelbrot_gpu(left, right, top, bottom);
	rect.setSize((sf::Vector2f)window->getSize());
	//rect.setTexture(&mandelbrot_tex);
}

// Render level
void Level::render()
{
	beginDraw();
	window->draw(rect);
	endDraw();
}

void Level::beginDraw()
{
	window->clear(sf::Color(100, 149, 237));
}

// Ends rendering to the back buffer, and swaps buffer to the screen.
void Level::endDraw()
{
	window->display();
}

void Level::compute_mandelbrot(double left, double right, double top, double bottom)
{
	sf::Vector2u Wsize = window->getSize();
	
	int colorPeriod = 0xFFFFFF / MAX_ITERATIONS;

	std::mutex pixelmute;
	
	for(int y = 0; y < Wsize.y; ++y)
	{
		std::cout << "Line " + std::to_string(y) + "\n";

		for (int x = 0; x < Wsize.x; ++x)
		{
			
			// Work out the point in the complex plane that
			// corresponds to this pixel in the output image.
			complex<double> c(left + (x * (right - left) / Wsize.x),
				top + (y * (bottom - top) / Wsize.y));

			// Start off z at (0, 0).
			complex<double> z(0.0, 0.0);

			// Iterate z = z^2 + c until z moves more than 2 units
			// away from (0, 0), or we've iterated too many times.
			sf::Uint32 iterations = 0;
			while (abs(z) < 2.0 && iterations < MAX_ITERATIONS)
			{
				z = (z * z) + c;

				++iterations;
			}
			//int r = (iterations > 510 ) ? iterations : 0;
			//int g = (iterations > 255 && iterations < 510) ? iterations : 0;
			//int b = (iterations < 255 ) ? iterations : 0;
			sf::Uint32 color = (0xFFFFFF-(iterations*colorPeriod)) | 255;
			
			mandelbrot->setPixel(x, y, sf::Color(color));

			//if (iterations == MAX_ITERATIONS)
			//{
	
			//	// z didn't escape from the circle.
			//	// This point is in the Mandelbrot set.
			//	mandelbrot.setPixel(x,y,sf::Color::Black); // black
			//}
			//else
			//{
			//	// z escaped within less than MAX_ITERATIONS
			//	// iterations. This point isn't in the set.
			//	mandelbrot.setPixel(x,y,sf::Color::White); // white
			//}
		}
	}
	mandelbrot_tex.update(*mandelbrot);
}

void Level::compute_mandelbrot_multicore(double left, double right, double top, double bottom)
{
	sf::Vector2u Wsize = window->getSize();
	
	
	int colorPeriod = 0xFFFFFF / MAX_ITERATIONS;

	std::vector<std::thread*> threads;

	//compute slice lambda
	auto compute_slice = [=](int start, int end)
	{
		for (int y = start; y < end; ++y)
		{
			std::cout << "Line " + std::to_string(y) + "\n";

			for (int x = 0; x < Wsize.x; ++x)
			{

				// Work out the point in the complex plane that
				// corresponds to this pixel in the output image.
				complex<double> c(left + (x * (right - left) / Wsize.x),
					top + (y * (bottom - top) / Wsize.y));

				// Start off z at (0, 0).
				complex<double> z(0.0, 0.0);

				// Iterate z = z^2 + c until z moves more than 2 units
				// away from (0, 0), or we've iterated too many times.
				sf::Uint32 iterations = 0;
				while (abs(z) < 2.0 && iterations < MAX_ITERATIONS)
				{
					z = (z * z) + c;

					++iterations;
				}
				//int r = (iterations > 510 ) ? iterations : 0;
				//int g = (iterations > 255 && iterations < 510) ? iterations : 0;
				//int b = (iterations < 255 ) ? iterations : 0;
				sf::Uint32 color = (0xFFFFFF - (iterations * colorPeriod)) | 255;


				mandelbrot->setPixel(x, y, sf::Color(color));

				//if (iterations == MAX_ITERATIONS)
				//{

				//	// z didn't escape from the circle.
				//	// This point is in the Mandelbrot set.
				//	mandelbrot.setPixel(x,y,sf::Color::Black); // black
				//}
				//else
				//{
				//	// z escaped within less than MAX_ITERATIONS
				//	// iterations. This point isn't in the set.
				//	mandelbrot.setPixel(x,y,sf::Color::White); // white
				//}
			}
		}
	};

	//threading
	int maxT = std::thread::hardware_concurrency();
	//cut into slices

	int slice = Wsize.y / maxT;

	//thread slices
	for (int i = 0; i < maxT; i++)
	{
		threads.push_back(new std::thread(compute_slice,i*slice,(i+1)*slice));
	}

	for (auto i : threads)
	{
		i->join();
	}
	
	//mandelbrot_tex = new sf::Texture;
	mandelbrot_tex.loadFromImage(*mandelbrot);

	delete mandelbrot;
	mandelbrot = NULL;
 }

void Level::compute_mandelbrot_gpu(double left, double right, double top, double bottom)
{
	const sf::Vector2i Wsize = { 1024,1024};
	int colorPeriod = 0xFFFFFF / MAX_ITERATIONS;

	//create array the size of window
	uint32_t* pImage = new uint32_t[Wsize.x * Wsize.y];
	//initialise to 0
	for (int i = 0; i < Wsize.x * Wsize.y; i++) pImage[i]=0;

	concurrency::extent<2> e (Wsize.x , Wsize.y);
	concurrency::array_view<uint32_t,2> a (e, pImage);

	a.discard_data();

	concurrency::parallel_for_each(a.extent, [=](concurrency::index<2> idx)restrict(amp)
	{
		//in 1D array size maxX*maxY any elements [index] x is index%maxX and y is index/maxX
		Complex1 c{ left + ((idx[0]) * (right - left) / Wsize.x),
			top + ((idx[1]) * (bottom - top) / Wsize.y)};
		Complex1 z{ 0.0,0.0 };

		uint32_t iterations = 0;

		while (c_abs(z) < 2.0 && iterations < MAX_ITERATIONS)
		{
			z = c_add(c_mul(z,z),c);
			++iterations;
		}
		
		uint32_t R;
		uint32_t G;
		uint32_t B;

		uint32_t color = ((iterations * colorPeriod)) | 255;
		a [idx] = color;
	});

	a.synchronize();

	for (int i = 0; i < Wsize.x; i++)
	{
		for (int j = 0; j < Wsize.y; j++)
		{
			if (mandelbrot)mandelbrot->setPixel(i, j, sf::Color((0xFFFFFF - pImage[(i * Wsize.y) + j])|255));

		}
	}

	//mandelbrot_tex.update(*mandelbrot);
	mandelbrot_tex.update(blur(mandelbrot));
	

	delete[] pImage;
	
}

//Image size must be 1024x1024
sf::Image Level::blur(sf::Image* image, int passes )
{
	if (passes <= 0) return *image;
	sf::Vector2i SIZE = { (int)image->getSize().x,(int)image->getSize().y };

	////Tilesizelimit 1024 total between xyz
	////no way to fit larger images without splitting them up or somehow overlaping the tiles
	//// 
	//// to blur an image this way would mean that the maximum resolution of an image would be 1024 x 1024 right?
	//// 
	////remove the tiles entirely or find a solution around that

	const int Filter_S = 7;

	uint32_t* iImage = new uint32_t[SIZE.x * SIZE.y];
	uint32_t* bImage = new uint32_t[SIZE.x * SIZE.y];

	for (int i = 0; i < SIZE.x; i++)
	{
		for (int j = 0; j < SIZE.y; j++)
		{
			iImage[(i * SIZE.y) + j] = image->getPixel(i, j).toInteger();
		}
	}

	static const int TSH = 1 << 10; //1024 maximum tile size
	static const int TSV = 1 << 10;

	//full filter
	struct filterA
	{
		float Distribution[7 * 7] =
		{
		  0.000904706, 0.003157733, 0.00668492,  0.008583607, 0.00668492,  0.003157733, 0.000904706,
		  0.003157733, 0.01102157,  0.023332663, 0.029959733, 0.023332663, 0.01102157,  0.003157733,
		  0.00668492,  0.023332663, 0.049395249, 0.063424755, 0.049395249, 0.023332663, 0.00668492,
		  0.008583607, 0.029959733, 0.063424755, 0.081438997, 0.063424755, 0.029959733, 0.008583607,
		  0.00668492,  0.023332663, 0.049395249, 0.063424755, 0.049395249, 0.023332663, 0.00668492,
		  0.003157733, 0.01102157,  0.023332663, 0.029959733, 0.023332663, 0.01102157,  0.003157733,
		  0.000904706, 0.003157733, 0.00668492,  0.008583607, 0.00668492,  0.003157733, 0.000904706
		};

	};

	//half matrix
	struct filterB {
		float Distribution[7] = { 1, 6 , 15, 20, 15, 6, 1 };
		float FRACTION = 1.f / 64; //1.f / 256.f;

		filterB() {
			for (auto i = 0; i < 7; i++)
			{
				Distribution[i] *= FRACTION;
			}
		}
	};

	filterB GD;
	extent<2> E(1024, 1024);

	concurrency::array_view<uint32_t, 2> av_src(E, iImage);
	concurrency::array_view<uint32_t, 2> av_dst(E, bImage);
	av_dst.discard_data();

	/*It is wise to use exception handling here - AMP can fail for many reasons
	and it useful to know why (e.g. using double precision when there is limited or no support)*/
	for (int i = 0; i < passes; i++) {
		try
		{
			concurrency::parallel_for_each(E.tile<TSH, 1>(), [=](tiled_index<TSH, 1> tidx) restrict(amp)
			{
				index<2> idx = tidx.global;

				tile_static uint32_t tile_data[TSH];
				tile_data[idx[0]] = av_src[idx];
				tidx.barrier.wait();

				int textureLocationX = idx[0] - (Filter_S / 2);

				uint32_t R = 0;
				uint32_t G = 0;
				uint32_t B = 0;

				for (int i = 0; i < Filter_S; i++)
				{
					R += ((tile_data[textureLocationX + i] & 0xFF000000) >> 24) * GD.Distribution[i];
					G += ((tile_data[textureLocationX + i] & 0x00FF0000) >> 16) * GD.Distribution[i];
					B += ((tile_data[textureLocationX + i] & 0x0000FF00) >> 8) * GD.Distribution[i];
				}
				if (R > 255) R = 255;
				if (G > 255) G = 255;
				if (B > 255) B = 255;

				tidx.barrier.wait();
				av_dst[idx] = (R << 24 | G << 16 | B << 8 | 255);

			});

			av_dst.synchronize();

			std::swap(av_src, av_dst);
			av_dst.discard_data();


			concurrency::parallel_for_each(E.tile<TSV, 1>(), [=](tiled_index<TSV, 1> tidx) restrict(amp)
			{
				index<2> idx = tidx.global;

				tile_static uint32_t tile_data[TSV];
				tile_data[idx[0]] = av_src[idx[1]][idx[0]];
				tidx.barrier.wait();

				int textureLocationX = idx[1] - (Filter_S / 2);

				uint32_t R = 0;
				uint32_t G = 0;
				uint32_t B = 0;

				for (int i = 0; i < Filter_S; i++)
				{
					R += ((tile_data[textureLocationX + i] & 0xFF000000) >> 24) * GD.Distribution[i];
					G += ((tile_data[textureLocationX + i] & 0x00FF0000) >> 16) * GD.Distribution[i];
					B += ((tile_data[textureLocationX + i] & 0x0000FF00) >> 8) * GD.Distribution[i];
				}

				if (R > 255) R = 255;
				if (G > 255) G = 255;
				if (B > 255) B = 255;

				tidx.barrier.wait();
				av_dst[idx[1]][idx[0]] = (R << 24 | G << 16 | B << 8 | 255);
			});

			av_dst.synchronize();
		}
		catch (const Concurrency::runtime_exception& ex)
		{
			MessageBoxA(NULL, ex.what(), "Error", MB_ICONERROR);
		}
	}

	sf::Image blurbedImage;

	blurbedImage.create(SIZE.x, SIZE.y, sf::Color::Cyan);

	for (int i = 0; i < SIZE.x; i++)
	{
		for (int j = 0; j < SIZE.y; j++)
		{
			blurbedImage.setPixel(i, j, sf::Color(bImage[(i * SIZE.y) + j]));
		}
	}

	delete[] bImage;
	delete[] iImage;

	return blurbedImage;
}



