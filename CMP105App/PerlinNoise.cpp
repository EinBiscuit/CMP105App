#include "PerlinNoise.h"

using namespace concurrency;

double* PerlinNoise::GeneratePerlin()
{
	PerlinClock::time_point start_0 = PerlinClock::now();

	srand(time(NULL));

	Seed seed_arr;

	double* result;

	result = new double[WIDTH * HEIGHT];
	for (int i = 0; i < WIDTH * HEIGHT; i++) result[i] = 0;

	extent<2> E(WIDTH, HEIGHT);
	array_view<double, 2> result_av(E, result);
	result_av.discard_data();

	//how many pixels are in a grid cell or essentially tile size
	//defined next to width and height in PerlinNoise.h

	PerlinClock::time_point start_1 = PerlinClock::now();

	try
	{
		concurrency::parallel_for_each(E, [=](index<2> idx) restrict(amp)
			{
				double fraction = 1. / PTS;

				//counterclockwise
				Vector2 GridCorners[4] = {

				seed_arr.Arr[idx[0] / PTS]       [idx[1] / PTS], //top left
				seed_arr.Arr[idx[0] / PTS]       [(idx[1] / PTS) + 1], //bottom left
				seed_arr.Arr[(idx[0] / PTS) + 1] [(idx[1] / PTS) + 1], //bottom right
				seed_arr.Arr[(idx[0] / PTS) + 1] [idx[1] / PTS] //top right

				//random number generation within thread a function
				//random(idx[0] / PTS,idx[1] / PTS),
				//random(idx[0] / PTS,(idx[1] / PTS) + 1),
				//random((idx[0] / PTS) + 1,(idx[1] / PTS) + 1),
				//random((idx[0] / PTS) + 1,idx[1] / PTS)
				};

				// displacement of the candidate point to the tile origin
				Vector2 point =
				{
					(double)(idx[0]%PTS)*fraction, (double)(idx[1]%PTS)*fraction
				};
				
				double dotone = dot(point , GridCorners[0]);
				double dottwo = dot({ point.x, point.y-1. }, GridCorners[1]);

				//__dp_d3d_smoothstepf produces a different result
				double smerp_0 = smoothstep(dotone ,dottwo, point.y);

				//do the x
				dotone = dot( { point.x-1., point.y-1. },GridCorners[2] );
				dottwo = dot({ point.x - 1., point.y }, GridCorners[3] );

				double smerp_1 = smoothstep(dotone, dottwo, point.y);

				result_av[idx] = smoothstep(smerp_0, smerp_1, point.x);
			});
	}
	catch (const Concurrency::runtime_exception& ex)
	{
		MessageBoxA(NULL, ex.what(), "Error", MB_ICONERROR);
	}

	result_av.synchronize();

	PerlinClock::time_point end = PerlinClock::now();
	// Compute the difference between the two times in milliseconds
	auto time = duration_cast<milliseconds>(end - start_1).count();
	std::cout << "Perlin Gpu Compute took: " << time << " ms." << std::endl;

	time = duration_cast<milliseconds>(end - start_0).count();
	std::cout << "Perlin Function took: " << time << " ms." << std::endl;

	//hard core debug
	/*for (int i = 0; i < HEIGHT; i++)
	{
		std::cout << "\nNEW ROW: ";
		for (int j = 0; j <WIDTH; j++)
		{
			std::cout <<" "<< result[i* WIDTH +j]<<" ";

		}

	}*/

	return result;
}

double* PerlinNoise::GeneratePerlinSerial()
{
	PerlinClock::time_point start = PerlinClock::now();

	srand(time(NULL));

	Seed seed_arr;

	double* result;

	result = new double[WIDTH * HEIGHT];
	for (int i = 0; i < WIDTH * HEIGHT; i++) result[i] = 0;

	//how many pixels are in a grid cell

	for(int x = 0 ; x < WIDTH; x++)
	{ 
		for(int y = 0; y < HEIGHT; y++)
		{
			double fraction = 1. / PTS;

			//counterclockwise
			Vector2 GridCorners[4] = {
				seed_arr.Arr[x / PTS]       [y / PTS], //top left
				seed_arr.Arr[x / PTS]       [y / PTS + 1], //bottom left
				seed_arr.Arr[x / PTS + 1] [y / PTS + 1], //bottom right
				seed_arr.Arr[x / PTS + 1] [y / PTS] //top right
			};


			// displacement of the candidate point to the tile origin
			Vector2 point =
			{
				(double)(x%PTS) * fraction, (double)(y%PTS) * fraction
			};

			double dotone = dot(point, GridCorners[0]);
			double dottwo = dot({ point.x, point.y - 1. }, GridCorners[1]);

			//__dp_d3d_smoothstepf produces a different result
			double smerp_0 = lerp(dotone, dottwo, point.y);

			//do the x
			dotone = dot({ point.x - 1., point.y - 1. }, GridCorners[2]);
			dottwo = dot({ point.x - 1., point.y }, GridCorners[3]);

			double smerp_1 = lerp(dotone, dottwo, point.y);

			result[(x*WIDTH)+y] = lerp(smerp_0, smerp_1, point.x);
		}
	}

	PerlinClock::time_point end = PerlinClock::now();
	// Compute the difference between the two times in milliseconds
	auto time = duration_cast<milliseconds>(end - start).count();
	std::cout << "Serial Perlin Function took: " << time << " ms." << std::endl;

	return result;
}

double* PerlinNoise::GeneratePerlin_cpu_fromWiki()
{
	PerlinClock::time_point start = PerlinClock::now();

	double* result;

	//to reproduce tiles
	int TS = 1 << 8;

	result = new double[WIDTH * HEIGHT];
	for (int i = 0; i < WIDTH * HEIGHT; i++) result[i] = 0;

	for (int x = 0; x < HEIGHT; x++) {
		for (int y = 0; y < WIDTH; y++) {

			int x0 = (x / TS);
			int x1 = x0 + 1;
			int y0 = (y / TS);
			int y1 = y0 + 1;

			// Determine interpolation weights
			// Could also use higher order polynomial/s-curve here
			float sx = ((float)x / TS) - x0;
			float sy = ((float)x / TS) - y0;

			//sx /= TS;
			//sy /= TS;

			// Interpolate between grid point gradients
			float n0, n1, ix0, ix1, value;

			n0 = dotGridGradient(x0, y0, x, y);
			n1 = dotGridGradient(x1, y0, x, y);
			ix0 = smoothstep(n0, n1, sx);

			n0 = dotGridGradient(x0, y1, x, y);
			n1 = dotGridGradient(x1, y1, x, y);
			ix1 = smoothstep(n0, n1, sx);

			result[x * WIDTH + y] = smoothstep(ix0, ix1, sy); // results to -1  1 range
		}

	}
	PerlinClock::time_point end = PerlinClock::now();
	auto time = duration_cast<milliseconds>(end - start).count();
	std::cout << "Serial Perlin Function took: " << time << " ms." << std::endl;

	return result;
}

