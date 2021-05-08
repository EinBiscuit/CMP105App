#pragma once
#include <amp.h>
#include <amp_math.h>
#include <random>
#include <math.h>
#include <iostream>
#include <chrono>

using std::chrono::duration_cast;
using std::chrono::milliseconds;

typedef std::chrono::steady_clock PerlinClock;

class PerlinNoise
{
public:
	// amp friendly vector 2 struct with 
	struct Vector2
	{
		double x;
		double y;
	};
	//FROM WIKIPEDIA  ADAPTED by me TO FOR AMP
	static Vector2 random(int x, int y) restrict(amp)
	{
		float random = 2920.f * concurrency::fast_math::sin(x * 21942.f + y * 171324.f + 8912.f) * concurrency::fast_math::cos(x * 23157.f * y * 217832.f + 9758.f);
		return { concurrency::fast_math::sin(random), concurrency::fast_math::cos(random) };
	}

	//amp friendly 
	static double dot(Vector2 a, Vector2 b) restrict(cpu, amp) {
		return ((a.x * b.x) + (a.y * b.y));
	}

	static double lerp(double a, double b, double w) restrict(cpu, amp) {
		return (b - a) * w + a;
	}

	static double smoothstep(double a, double b, double w) restrict(cpu, amp) {
		return (b - a) * (3.0 - w * 2.0) * w * w + a;
	}

	static void normalise(Vector2& a)
	{
		float v = (a.x * a.x) + (a.y * a.y);
		if (v > 0) {
			v = sqrtf(v);
			a.x /= v;
			a.y /= v;
		};
	}

	// AMD SMOOTHSTEP
	static float smoothstepAMD(float edge0, float edge1, float x) restrict(cpu, amp)
	{
		// Scale, bias and saturate x to 0..1 range
		x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
		// Evaluate polynomial
		return x * x * (3 - 2 * x);
	}

	static float clamp(float x, float lowerlimit, float upperlimit) restrict(cpu, amp)
	{
		if (x < lowerlimit)
			x = lowerlimit;
		if (x > upperlimit)
			x = upperlimit;
		return x;
	}


	// FROM WIKIPEDIA https://en.wikipedia.org/wiki/Perlin_noise
	static Vector2 randomGradient(int ix, int iy) {
		// Random float. No precomputed gradients mean this works for any number of grid coordinates
		float random = 2920.f * sin(ix * 21942.f + iy * 171324.f + 8912.f) * cos(ix * 23157.f * iy * 217832.f + 9758.f);
		return { cos(random), sin(random) };
	}

	static float dotGridGradient(int ix, int iy, float x, float y) {
		// Get gradient from integer coordinates
		Vector2 gradient = randomGradient(ix, iy);

		normalise(gradient);

		// Compute the distance vector
		float dx = x - (float)ix;
		float dy = y - (float)iy;

		// Compute the dot-product
		return (dx * gradient.x + dy * gradient.y);
	}
	//

	//amp friendly seed vector array
	struct Seed
	{
		Vector2 Arr[16][16];
		Seed() {
			int x = 16;
			int y = 16;
			for (int i = 0; i < x; i++) {
				for (int j = 0; j < y; j++) {

					int xneg = rand() % 2;
					int yneg = rand() % 2;

					Arr[i][j] = { (xneg) ? -(double)rand() : (double)rand(), (yneg) ? -(double)rand() : (double)rand() };
					normalise(Arr[i][j]);
				}
			}
		}
	};

	static const uint32_t WIDTH = 1 << 10;
	static const uint32_t HEIGHT = 1 << 10;

	static double* GeneratePerlin();
	//very noisy, work but i think it mostly displays the random function
	static double* GeneratePerlin_cpu_fromWiki();

	//void GeneratePerlinLite();
};

