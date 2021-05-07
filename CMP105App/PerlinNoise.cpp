#include "PerlinNoise.h"

using namespace concurrency;

double* PerlinNoise::GeneratePerlin()
{
	srand(time(NULL));

	Seed seed_arr;

	double* result;

	result = new double[WIDTH * HEIGHT];
	for (int i = 0; i < WIDTH * HEIGHT; i++) result[i] = 0;

	extent<2> E(WIDTH, HEIGHT);
	array_view<double, 2> result_av(E, result);
	result_av.discard_data();

	//how many pixels are in a grid cell
	const static int PTS = 64;
	try
	{
		concurrency::parallel_for_each(E.tile<PTS, 2>(), [=](tiled_index<PTS, 2> tidx) restrict(amp)
			{
				index<2> idx = tidx.global;

				//counterclockwise
				Vector2 GridCorners[4] = {

				seed_arr.Arr[(int)idx[0] / PTS][idx[1] / PTS],
				seed_arr.Arr[(int)idx[0] / PTS][(idx[1] / PTS) + 1],
				seed_arr.Arr[((int)idx[0] / PTS) + 1][(idx[1] / PTS) + 1],
				seed_arr.Arr[((int)idx[0] / PTS) + 1][idx[1] / PTS]

				//random(idx[0] / PTS,idx[1] / PTS),
				//random(idx[0] / PTS,(idx[1] / PTS) + 1),
				//random((idx[0] / PTS) + 1,(idx[1] / PTS) + 1),
				//random((idx[0] / PTS) + 1,idx[1] / PTS)
				};



				Vector2 point = { (double)tidx.local[0],(double)tidx.local[1] };
				//normalize the vector
				double squaresum = (point.x * point.x) + (point.y * point.y);
				if (squaresum > 0)
				{
					squaresum = fast_math::sqrt(squaresum);
					point = { point.x / squaresum , point.y / squaresum };
				}
				//do the y
				double dotone = dot(GridCorners[0], point);
				double dottwo = dot(GridCorners[1], { point.x,-point.y });

				//__dp_d3d_smoothstepf produces a different result
				double smerp_y = lerp(dotone, dottwo, (float)tidx.local[1] / PTS);

				//do the x
				dotone = dot(GridCorners[2], { -point.x,-point.y });
				dottwo = dot(GridCorners[3], { -point.x, point.y });

				double smerp_x = lerp(dotone, dottwo, (float)tidx.local[1] / PTS);

				result_av[idx] = lerp(smerp_y, smerp_x, (float)tidx.local[0] / PTS);
			});
	}
	catch (const Concurrency::runtime_exception& ex)
	{
		MessageBoxA(NULL, ex.what(), "Error", MB_ICONERROR);
	}

	result_av.synchronize();

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

double* PerlinNoise::GeneratePerlin_cpu_fromWiki()
{
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
			ix0 = lerp(n0, n1, sx);

			n0 = dotGridGradient(x0, y1, x, y);
			n1 = dotGridGradient(x1, y1, x, y);
			ix1 = lerp(n0, n1, sx);

			result[x * WIDTH + y] = lerp(ix0, ix1, sy); // results to -1  1 range
		}

	}

	return result;
}
