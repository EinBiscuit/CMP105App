#pragma once
#include <amp.h>
#include <amp_math.h>
#include <random>
#include <math.h>
#include <iostream>

class PerlinNoise
{
public:
	// amp friendly vector 2 struct with 
	struct Vector2
	{
		double x;
		double y;
	};

	static double dot(Vector2 a, Vector2 b) restrict(cpu, amp) {
		return ((a.x * b.x) + (a.y* b.y));
	}

	static double lerp(double a, double b, double w) restrict(cpu, amp) {
		return (b - a) * w + a;
	}

	static double smoothstep(double a, double b, double w) restrict(cpu, amp) {
		return (b - a) * (3.0 - w * 2.0) * w * w + a;
	}

	static void normalise(Vector2 a)
	{
		double v = sqrtf(a.x * a.x + a.y * a.y);
		a.x /= v;
		a.y /= v;
	}

	/*static Vector2 lerp(Vector2 a, Vector2 b, double factor) restrict(cpu, amp) 
	{
		return { a + ((b - a) * factor)
	};
	*/

	//amp friendly seed vector array
	struct Seed
	{
		Vector2 Arr[16][16];
		Seed() {
			int x = 16;
			int y = 16;
			for (int i = 0; i < x; i++) {
				for (int j = 0; j < y; j++) {
					Arr[i][j] = { (float)(rand() % 100), (float)(rand() % 100) };
					normalise(Arr[i][j]);
				}
			}
		}
	};

	static const uint32_t WIDTH = 1 << 8;
	static const uint32_t HEIGHT = 1 << 8;

	static double* GeneratePerlin();
	//void GeneratePerlinLite();
};

