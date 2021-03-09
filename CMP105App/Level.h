#pragma once

#include <SFML/Graphics.hpp>
#include <string.h>
#include <iostream>


struct Mouse;

class Level{
public:
	Level(sf::RenderWindow* hwnd, Mouse* mus);
	~Level();
	void compute_mandelbrot(double left, double right, double top, double bottom);
	void compute_mandelbrot_multicore(double left, double right, double top, double bottom);
	void compute_mandelbrot_gpu(double left, double right, double top, double bottom);
	void handleInput();
	void update();
	void render();

private:
	void beginDraw();
	void endDraw();
	sf::RenderWindow* window;
	Mouse* mouse;

	sf::RectangleShape rect;
	sf::Texture mandelbrot_tex;

	int depth = 1;
	sf::Vector2f center;
	double left = -2.0;
	double right = 1.0;
	double top = 1.125;
	double bottom = -1.125;

	
};