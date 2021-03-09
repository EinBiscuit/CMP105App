#include "Level.h"
#include <math.h>
#include <cstdint>
#include <cstdlib>
#include <complex>
#include <thread>
#include <mutex>
#include <amp.h>
#include <amp_math.h>

#define MAX_ITERATIONS 80

using namespace concurrency;
using std::complex;


struct Complex1 { float x; float y; };

Complex1 c_add(Complex1 c1, Complex1 c2) restrict(cpu, amp) // restrict keyword -able to execute this function on the GPU and CPU
{
	Complex1 tmp;
	float a = c1.x;
	float b = c1.y;
	float c = c2.x;
	float d = c2.y;
	tmp.x = a + c;
	tmp.y = b + d;
	return tmp;
}
//c_add
float c_abs(Complex1 c) restrict(cpu, amp)
{
	return concurrency::fast_math::sqrt(c.x * c.x + c.y * c.y);
}
// c_abs 
Complex1 c_mul(Complex1 c1, Complex1 c2) restrict(cpu, amp)
{
	Complex1 tmp;
	float a = c1.x;
	float b = c1.y;
	float c = c2.x;
	float d = c2.y;
	tmp.x = a * c - b * d;
	tmp.y = b * c + a * d;
	return tmp;
} // c_mul

Level::Level(sf::RenderWindow* hwnd, Mouse* mus)
{
	window = hwnd;
	mouse = mus;
	// initialise game objects

	// instead of running entire queery_amp i just ran first line of it
	std::cout << accelerator::get_all().size() << std::endl;

	//compute_mandelbrot(-0.751085, -0.734975, 0.118378, 0.134488);
	//compute_mandelbrot_multicore(-2.0, 1.0, 1.125, -1.125);
	compute_mandelbrot_gpu(-2.0, 1.0, 1.125, -1.125);
	
	//rectangle
	rect.setSize((sf::Vector2f)window->getSize());
	rect.setTexture(&mandelbrot_tex);
	//rect.setFillColor(sf::Color::Red);
}

Level::~Level()
{
}

// handle user input
void Level::handleInput()
{
  
}

// Update game objects
void Level::update()
{


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
	sf::Image mandelbrot;
	mandelbrot.create(Wsize.x, Wsize.y, sf::Color(0));
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
			
			mandelbrot.setPixel(x, y, sf::Color(color));

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

	mandelbrot_tex.loadFromImage(mandelbrot);
}

void Level::compute_mandelbrot_multicore(double left, double right, double top, double bottom)
{
	sf::Vector2u Wsize = window->getSize();
	sf::Image* mandelbrot;
	mandelbrot = new sf::Image;
	mandelbrot->create(Wsize.x, Wsize.y, sf::Color(0));
	
	int colorPeriod = 0xFFFFFF / MAX_ITERATIONS;

	std::vector<std::thread*> threads;

	//compute slice lambda
	auto compute_slice = [Wsize,left,right,top,bottom,colorPeriod, mandelbrot](int start, int end)
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
		
	mandelbrot_tex.loadFromImage(*mandelbrot);

	delete mandelbrot;
	mandelbrot = NULL;
 }

void Level::compute_mandelbrot_gpu(double left, double right, double top, double bottom)
{
	const sf::Vector2u Wsize = window->getSize();
	int colorPeriod = 0xFFFFFF / MAX_ITERATIONS;

	uint32_t* pImage = new uint32_t[(uint64_t)Wsize.x * (uint64_t)Wsize.y];
	for (int i = 0; i < Wsize.x * Wsize.y; i++) pImage[i] = 0;

	concurrency::extent<1> e (Wsize.x * Wsize.y);
	concurrency::array_view<uint32_t> a (e, pImage);

	a.discard_data();

	concurrency::parallel_for_each(a.extent, [=](concurrency::index<1> idx)restrict(amp)
		{
			Complex1 c{ left + ((idx[0]%Wsize.y) * (right - left) / Wsize.x),
				top + ((idx[0]/Wsize.x) * (bottom - top) / Wsize.y) };
			Complex1 z{ 0.0,0.0 };

			uint32_t iterations = 0;

			while (c_abs(z) < 2.0 && iterations < MAX_ITERATIONS)
			{
				z = c_add(c_mul(z,z),c);
				++iterations;
			}

			uint32_t color = (0xFFFFFF - (iterations * colorPeriod)) | 255;

			a[idx] = color;
		});

	a.synchronize();

	sf::Image* mandelbrot;

	mandelbrot = new sf::Image;
	mandelbrot->create(Wsize.x,Wsize.y,sf::Color(0));

	for (int i = 0; i < Wsize.x; i++)
	{
		for (int j = 0; j < Wsize.y; j++)
		{
			mandelbrot->setPixel(i,j,sf::Color(pImage[i*Wsize.x + j]));
		}
	}

	mandelbrot_tex.loadFromImage(*mandelbrot);

	delete mandelbrot;
	mandelbrot = NULL;
}