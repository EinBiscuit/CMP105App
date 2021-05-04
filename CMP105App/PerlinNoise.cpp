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

	concurrency::parallel_for_each(E.tile<PTS, 2>(), [=](tiled_index<PTS, 2> tidx) restrict(amp)
		{
			index<2> idx = tidx.global;

			//counterclockwise
			Vector2 GridCorners[4] = {
			seed_arr.Arr[idx[0] / PTS][idx[1] / PTS],
			seed_arr.Arr[idx[0] / PTS][(idx[1] / PTS) + 1],
			seed_arr.Arr[(idx[0] / PTS) + 1][(idx[1] / PTS) + 1],
			seed_arr.Arr[(idx[0] / PTS) + 1][idx[1] / PTS]
			};

			Vector2 point = { (double)tidx.local[0],(double)tidx.local[1] };
			//normalize the vector 
			double len = fast_math::sqrtf(point.x * point.x + point.y * point.y);
			point = {point.x/len,point.y/len};
			//do the y
			double dotone = dot(GridCorners[0], point);
			double dottwo = dot(GridCorners[1], { point.x,-point.y });

			//__dp_d3d_smoothstepf produces a different result
			double smerp_y = smoothstep(dotone, dottwo, (double)tidx.local[1]/(double)PTS);

			//do the x
			dotone = dot(GridCorners[2], { -point.x,-point.y });
			dottwo = dot(GridCorners[3], { -point.x,point.y });

			
			float smerp_x = smoothstep(dotone, dottwo, (double)tidx.local[1]/(double)PTS);

			result_av[idx] = smoothstep(smerp_y, smerp_x, (double)tidx.local[0]/(double)PTS);
		});

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
