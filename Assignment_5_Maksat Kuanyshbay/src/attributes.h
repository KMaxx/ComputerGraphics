#pragma once

#include <Eigen/Core>

struct Light {
	Eigen::Vector3d position;
	Eigen::Vector3d intensity;
};

struct Camera {
	bool is_perspective;
	Eigen::Vector3d position;
	double field_of_view; // between 0 and PI
	double focal_length;
	double lens_radius; // for depth of field
};

struct Material {
	Eigen::Vector3d ambient_color;
	Eigen::Vector3d diffuse_color;
	Eigen::Vector3d specular_color;
	double specular_exponent; // Also called "shininess"

	Eigen::Vector3d reflection_color;
	Eigen::Vector3d refraction_color;
	double refraction_index;
};

class VertexAttributes
{
	public:
	VertexAttributes(float x = 0, float y = 0, float z = 0, float w = 1)
	{
		position << x,y,z,w;
		color << 1.0,1.0,1.0,1.0;
		normal << x,y,z,w;
		num = 0.0;
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
        // r.position = alpha*(a.position/a.position[3]) + beta*(b.position/b.position[3]) + gamma*(c.position/c.position[3]);
		// r.color = alpha*(a.color/a.color[3]) + beta*(b.color/b.color[3]) + gamma*(c.color/c.color[3]);
		// r.normal = alpha*(a.normal/a.normal[3]) + beta*(b.normal/b.normal[3]) + gamma*(c.normal/c.normal[3]);
	    r.position = alpha*a.position + beta*b.position + gamma*c.position;
		r.color = alpha*a.color + beta*b.color + gamma*c.color;
		r.normal = alpha*a.normal + beta*b.normal + gamma*c.normal;
        
		
		
		
		return r;
    }

	Eigen::Vector4f position;
	Eigen::Vector4f color;
	Eigen::Vector4f normal;
	float num;
};

class FragmentAttributes
{
	public:
	FragmentAttributes(float r = 0, float g = 0, float b = 0, float a = 1)
	{
		color << r,g,b,a;
	}

	Eigen::Vector4f color;
	Eigen::Vector4f position;
};

class FrameBufferAttributes
{
	public:
	FrameBufferAttributes(uint8_t r = 0, uint8_t g = 0, uint8_t b = 0, uint8_t a = 255)
	{
		color << r,g,b,a;
		depth = 2; 
	}

	Eigen::Matrix<uint8_t,4,1> color;
	float depth;
};

class UniformAttributes
{
	public:
	Eigen::Matrix4f projection;
	Material material;
	Camera camera;
	Light light;
	Eigen::Vector4f color;
};