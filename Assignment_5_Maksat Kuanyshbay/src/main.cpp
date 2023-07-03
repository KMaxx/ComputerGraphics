// C++ include
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Eigen for matrix operations
#include <Eigen/Dense>
#include <Eigen/Geometry>

// Utilities for the Assignment
#include "raster.h"
#include <gif.h>
// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

using namespace std;

FrameBuffer scale_down_4x(const FrameBuffer& fb)
{
	// The size of the framebuffer must be a multiple of 4
	assert(fb.rows() % 4 == 0);
	assert(fb.cols() % 4 == 0);

	// Allocate the reduced buffer
	FrameBuffer out(fb.rows()/4,fb.cols()/4);

	for (unsigned i=0;i<out.rows();i++)
	{
		for (unsigned j=0;j<out.cols();j++)
		{
			Eigen::Vector4f avg = Eigen::Vector4f::Zero();
			for (unsigned ii=0;ii<4;ii++)
				for (unsigned jj=0;jj<4;jj++)
					avg += fb(i*4+ii,j*4+jj).color.cast<float>();
			avg /= 16;
			out(i,j).color = avg.cast<uint8_t>();
		}
	}
	return out;
}


int main() 
{

	// The Framebuffer storing the image rendered by the rasterizer
	Eigen::Matrix<FrameBufferAttributes,Eigen::Dynamic,Eigen::Dynamic> frameBuffer(500,500);

	// Global Constants (empty in this example)
	UniformAttributes uniform;

	// Basic rasterization program
	Program program;

	// The vertex shader is the identity
	program.VertexShader = [](const VertexAttributes& va, const UniformAttributes& uniform)
	{
		VertexAttributes out;
		out.position = uniform.projection * va.position;
		out.normal = uniform.projection * va.normal;
		// out.position = uniform.projection * va.position;
		// out.position = out.position/out.position(3);
		// out.normal = uniform.projection * va.normal;
		// out.normal = out.normal/out.normal(3);

		Eigen::Vector3f hitposition = va.position.head<3>();
		Eigen::Vector3f lightposition(uniform.light.position(0), uniform.light.position(1), uniform.light.position(2));
		Eigen::Vector3f Li = (lightposition - hitposition).normalized();
		Eigen::Vector3f N = va.normal.head<3>();
		Eigen::Vector3f diffuse(uniform.material.diffuse_color(0), uniform.material.diffuse_color(1), uniform.material.diffuse_color(2));
		diffuse = diffuse * max(2.0*Li.dot(N),0.0);
		Eigen::Vector3f Hi = (Li-Eigen::Vector3f(0,0,-1)).normalized();
		Eigen::Vector3f specular(uniform.material.specular_color(0), uniform.material.specular_color(1), uniform.material.specular_color(2));
		specular = specular * pow(max(N.dot(Hi), (float)0.0), (float)1.0);
		Eigen::Vector3f ambient_color(uniform.material.ambient_color(0),uniform.material.ambient_color(1),uniform.material.ambient_color(2));
		Eigen::Vector3f ambient_light(0.7,0.7,0.7);
		ambient_color = ambient_color.array() * ambient_light.array();
		Eigen::Vector3f color = ambient_color + (diffuse + specular).cwiseProduct(Eigen::Vector3f(0.06,0.06,0.06)) /  (lightposition - hitposition).squaredNorm();
		out.color = Eigen::Vector4f(color(0),color(1),color(2), (float)1);

		return out;
	};

	// The fragment shader uses a fixed color
	program.FragmentShader = [](const VertexAttributes& va, const UniformAttributes& uniform)
	{

		FragmentAttributes out(va.color[0],va.color[1],va.color[2],va.color[3]);
		out.position = va.position;
		return out;
	};

	// The blending shader converts colors between 0 and 1 to uint8
	program.BlendingShader = [](const FragmentAttributes& fa, const FrameBufferAttributes& previous)
	{
		if (fa.position[2] < previous.depth)
		{
			FrameBufferAttributes out(fa.color[0]*255, fa.color[1]*255, fa.color[2]*255, fa.color[3]*255);
			out.depth = fa.position[2];
			return out;
		}
		else
			return previous;
	};

	string filename = string(DATA_DIR) + string("bunny.off");

	ifstream in(filename);
	string token;
	in >> token;
	int nv, nf, ne;
	in >> nv >> nf >> ne;

	Eigen::MatrixXd V(nv, 3); // n x 3 matrix (n points)
	Eigen::MatrixXf Normals(nv,4);
	Normals.setZero();
	Eigen::MatrixXi F(nf, 3); // m x 3 matrix (m triangles)
	double l = numeric_limits<int>::max();
	double b = numeric_limits<int>::max();
	double n = numeric_limits<int>::max();
	double r = -numeric_limits<int>::max();
	double t = -numeric_limits<int>::max();
	double f = -numeric_limits<int>::max();
	for (int i = 0; i < nv; ++i) {
		in >> V(i, 0) >> V(i, 1) >> V(i, 2);
		l = min(l, V(i,0));
		b = min(b, V(i,1));
		n = min(n, V(i,2));
		r = max(r, V(i,0));
		t = max(t, V(i,1));
		f = max(f, V(i,2));
	}
	
	for (int i = 0; i < nf; ++i) {
		int s;
		in >> s >> F(i, 0) >> F(i, 1) >> F(i, 2);
		assert(s == 3);
	}

	// One triangle in the center of the screen
	vector<VertexAttributes> vertices;

	vector<VertexAttributes> wvertices;
	for (int i = 0; i < F.rows(); ++i) {
	 	for (int k = 0; k < F.cols(); ++k) {
	 		vertices.push_back(VertexAttributes(V.row(F(i, k))(0),V.row(F(i, k))(1),V.row(F(i, k))(2)));
	 	}
		wvertices.push_back(VertexAttributes(V.row(F(i, 0))(0),V.row(F(i, 0))(1),V.row(F(i, 0))(2),1.0));
		wvertices.push_back(VertexAttributes(V.row(F(i, 1))(0),V.row(F(i, 1))(1),V.row(F(i, 1))(2),1.0));
		wvertices.push_back(VertexAttributes(V.row(F(i, 1))(0),V.row(F(i, 1))(1),V.row(F(i, 1))(2),1.0));
		wvertices.push_back(VertexAttributes(V.row(F(i, 2))(0),V.row(F(i, 2))(1),V.row(F(i, 2))(2),1.0));
		wvertices.push_back(VertexAttributes(V.row(F(i, 2))(0),V.row(F(i, 2))(1),V.row(F(i, 2))(2),1.0));
		wvertices.push_back(VertexAttributes(V.row(F(i, 0))(0),V.row(F(i, 0))(1),V.row(F(i, 0))(2),1.0));

		vertices[vertices.size()-3].normal << (vertices[vertices.size()-2].position.head<3>()-vertices[vertices.size()-3].position.head<3>()).cross(vertices[vertices.size()-1].position.head<3>()-vertices[vertices.size()-3].position.head<3>()).normalized(),1;
		// vertices[vertices.size()-2].normal = vertices[vertices.size()-3].normal;
		// vertices[vertices.size()-1].normal = vertices[vertices.size()-3].normal;
 		Normals.row(F(i,0)).transpose() = Normals.row(F(i,0)).transpose() + vertices[vertices.size()-3].normal;
		Normals.row(F(i,1)).transpose() = Normals.row(F(i,1)).transpose() + vertices[vertices.size()-3].normal;
		Normals.row(F(i,2)).transpose() = Normals.row(F(i,2)).transpose() + vertices[vertices.size()-3].normal;
		vertices[vertices.size()-3].normal = Normals.row(F(i,0)).transpose();

		//vertices[vertices.size()-2].normal << (vertices[vertices.size()-2].position.head<3>()-vertices[vertices.size()-3].position.head<3>()).cross(vertices[vertices.size()-1].position.head<3>()-vertices[vertices.size()-3].position.head<3>()).normalized(),1;

		vertices[vertices.size()-2].normal = Normals.row(F(i,1)).transpose();		
	
		vertices[vertices.size()-1].normal = Normals.row(F(i,2)).transpose();	

		// wvertices[wvertices.size()-6].normal = vertices[vertices.size()-3].normal;
		// wvertices[wvertices.size()-5].normal = vertices[vertices.size()-3].normal;
		// wvertices[wvertices.size()-3].normal = vertices[vertices.size()-3].normal;
		// wvertices[wvertices.size()-1].normal = vertices[vertices.size()-3].normal;
		// wvertices[wvertices.size()-2].normal = vertices[vertices.size()-3].normal;
		// wvertices[wvertices.size()-4].normal = vertices[vertices.size()-3].normal;
	}
	// Add a projective transformation
	for (int i = 0; i < vertices.size(); ++i) {
		vertices[i].normal = vertices[i].normal.normalized();
	}

	uniform.light.position << 0.5,-0.2,f+n;
	uniform.light.intensity << 10,10,10;


	uniform.camera.is_perspective = false;
	//uniform.camera.position << (float)(l+r)/2,(float)(b+t)/2,n-0.1;
	uniform.camera.position << 0,0,(n+f);
	uniform.camera.field_of_view = 0.7854;
	uniform.camera.focal_length = 5.0;
	uniform.camera.lens_radius = 0.08;

	uniform.material.ambient_color << 0.0, 0.5, 0.0;
	uniform.material.diffuse_color << 0.5, 0.5, 0.5;
	uniform.material.specular_color << 0.2, 0.2, 0.2;
	uniform.material.specular_exponent = 1.0;
	uniform.material.reflection_color << 0.7, 0.7, 0.7;
	uniform.material.refraction_color << 0.0, 0.0, 0.0;
	uniform.material.refraction_index = 1.0;
	uniform.color << 0.0, 0.0, 0.0,0.0;

	Eigen::Matrix4f m_pers;
	Eigen::Matrix4f m_orth;
	// float aspect_ratio = float(frameBuffer.cols())/float(frameBuffer.rows());
	// t = tan(uniform.camera.field_of_view)*n;
	// r = t*aspect_ratio;
	// l = -r;
	// b = -t;



	m_pers <<
	n, 0.0, 0.0, 0.0,
	0.0, n, 0.0, 0.0,
	0.0, 0.0, n+f, (-1.0*f*n),
	0.0, 0.0, 1.0, 0.0;	
	
	m_orth <<
	2/(r-l), 0.0, 0.0, -(r+l)/(r-l),
	0.0, 2/(t-b), 0.0, -(t+b)/(t-b),
	0.0, 0.0, 2/(n-f), -(n+f)/(n-f),
	0.0, 0.0, 0.0, 1.0;	
	
	Eigen::Vector3d e = uniform.camera.position;
	Eigen::Vector3d g = Eigen::Vector3d(0,0,-1);
	Eigen::Vector3d tt = Eigen::Vector3d(0,1,0);
	Eigen::Vector3d w = -g.normalized();
	Eigen::Vector3d u = tt.cross(w).normalized();
	Eigen::Vector3d v = w.cross(u);

	Eigen::Matrix4f m_cam;
	m_cam <<
	u(0), v(0), w(0), e(0),
	u(1), v(1), w(1), e(1),
	u(2), v(2), w(2), e(2),
	0.0, 0.0, 0.0, 1.0;

	Eigen::Vector4f bc;
	bc<<(float)(l+r)/2,(float)(b+t)/2,(float)(n+f)/2,1.0;
	Eigen::Matrix4f m_tb;
	m_tb <<
	1.0, 0.0, 0.0, bc(0),
	0.0, 1.0, 0.0, bc(1),
	0.0, 0.0, 1.0, bc(2),
	0.0, 0.0, 0.0, 1.0;
	Eigen::Matrix4f m_tbb;
	m_tbb <<
	1.0, 0.0, 0.0, -bc(0),
	0.0, 1.0, 0.0, -bc(1),
	0.0, 0.0, 1.0, -bc(2),
	0.0, 0.0, 0.0, 1.0;



	const char * fileName = "triangle.gif";
	vector<uint8_t> image;
	int delay = 25;
	GifWriter gi;
	GifBegin(&gi, fileName, frameBuffer.rows(), frameBuffer.cols(), delay);

	for (float i=0;i<2*3.14;i+=0.2512)
	{
		frameBuffer.setConstant(FrameBufferAttributes());
		Eigen::Matrix4f m_tt;
		m_tt <<
		1.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, +i*0.05,
		0.0, 0.0, 0.0, 1.0;
		Eigen::Matrix4f m_rot;
		m_rot <<
		cos(i), 0, -sin(i), 0,
		0, 1, 0, 0,
		sin(i), 0, cos(i), 0,
		0, 0, 0, 1;	
		uniform.projection = m_orth*m_cam.inverse()*m_tt*m_tb*m_rot*m_tbb;
		//uniform.projection = m_orth*m_cam.inverse()*m_tb*m_rot*m_tbb;
		rasterize_triangles(program,uniform,vertices,frameBuffer);
		//rasterize_lines(program,uniform,wvertices,0.8,frameBuffer);
		framebuffer_to_uint8(frameBuffer,image);
		GifWriteFrame(&gi, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
	}

	GifEnd(&gi);

	// uniform.projection = m_orth*m_pers*m_cam.inverse();
	// uniform.projection = m_orth*m_cam.inverse();
	// rasterize_triangles(program,uniform,vertices,frameBuffer);
	// rasterize_lines(program,uniform,wvertices,0.8,frameBuffer);

	// vector<uint8_t> image;
	// //frameBuffer = scale_down_4x(frameBuffer);
	// framebuffer_to_uint8(frameBuffer,image);
	// stbi_write_png("triangle.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows()*4);
	
	return 0;
}
