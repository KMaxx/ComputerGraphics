#pragma once

#include <Eigen/Core>

class VertexAttributes
{
	public:
	VertexAttributes(float x = 0, float y = 0, float z = 0, float w = 1)
	{
		position << x,y,z,w;
		color << 1,1,1,1;
		index = -1;
	}

    // Interpolates the vertex attributes
    static VertexAttributes interpolate(
        const VertexAttributes& a,
        const VertexAttributes& b,
        const VertexAttributes& c,
        const float alpha, 
        const float beta, 
        const float gamma
    ) 
    {
        VertexAttributes r;
        r.position = alpha*a.position + beta*b.position + gamma*c.position;
		r.color = alpha*a.color + beta*b.color + gamma*c.color;
        return r;
    }

	Eigen::Vector4f position;
	Eigen::Vector4f color;
	int index;
};

class FragmentAttributes
{
	public:
	FragmentAttributes(float r = 0, float g = 0, float b = 0, float a = 1)
	{
		color << r,g,b,a;
	}

	Eigen::Vector4f color;
};

class FrameBufferAttributes
{
	public:
	FrameBufferAttributes(uint8_t r = 0, uint8_t g = 0, uint8_t b = 0, uint8_t a = 255)
	{
		color << r,g,b,a;
	}

	Eigen::Matrix<uint8_t,4,1> color;
};

class UniformAttributes
{
	public:
		int clicks;
		float width;
		float height;
		std::vector<Eigen::Matrix4f> trans;
		std::vector<Eigen::Matrix4f> keyfs;
		int a;
		int b;
		int c;
		int col;
		Eigen::Vector4f acolor;
		Eigen::Vector4f bcolor;
		Eigen::Vector4f ccolor;
		Eigen::Vector4f p;
		Eigen::Vector2f s;
		bool is_animation;
		bool is_bezier;
		bool is_pressed;
		char mode;
		Eigen::Matrix4f zoom;
		Eigen::Matrix4f view;
		Eigen::Vector4f bc;
		Eigen::Matrix4f tc;
		Eigen::Matrix4f ct;
		Eigen::Matrix4f rt;
};