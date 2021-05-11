#pragma once
#include <amp_math.h>


struct Complex1 { double x; double y; };

Complex1 c_add(Complex1 c1, Complex1 c2) restrict(cpu, amp) // restrict keyword -able to execute this function on the GPU and CPU
{
	Complex1 tmp;
	double a = c1.x;
	double b = c1.y;
	double c = c2.x;
	double d = c2.y;
	tmp.x = a + c;
	tmp.y = b + d;
	return tmp;
}
//c_add
double c_abs(Complex1 c) restrict(cpu, amp)
{
	return concurrency::fast_math::sqrt((float)(c.x * c.x + c.y * c.y));

}
// c_abs 
Complex1 c_mul(Complex1 c1, Complex1 c2) restrict(cpu, amp)
{
	Complex1 tmp;
	double a = c1.x;
	double b = c1.y;
	double c = c2.x;
	double d = c2.y;
	tmp.x = a * c - b * d;
	tmp.y = b * c + a * d;
	return tmp;
} // c_mul

