// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"


// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

Vector3d intersection_sphere(Vector3d e, Vector3d d, Vector3d c, double r) {
        double A = d.transpose()*d;
		double B = 2*d.transpose()*(e-c);
		double C = (e-c).transpose()*(e-c) - r*r;
		double disc = B*B-4*A*C;
		Vector3d x(0,0,-1);
		if(disc>0){
			x << (-B - sqrt(disc))/(2*A), (-B + sqrt(disc))/(2*A), 1;
		}
		return x;
    }

void raytrace_sphere() {
	std::cout << "Simple ray tracer, one sphere with orthographic projection" << std::endl;

	const std::string filename("sphere_orthographic.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// Single light source
	const Vector3d light_position(-1,1,1);

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// Prepare the ray
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = RowVector3d(0,0,-1);

			// Intersect with the sphere
			Vector2d ray_on_xy(ray_origin(0),ray_origin(1));
			const Vector3d sphere_origin(0,0,0);
			const double sphere_radius = 0.9;
			Vector3d x = intersection_sphere(ray_origin, ray_direction, sphere_origin, sphere_radius);

			if (x(2) > 0) {
				// The ray hit the sphere, compute the exact intersection point
				Vector3d ray_intersection = ray_origin + std::min(x(0),x(1))*ray_direction;

				// Compute normal at the intersection point
				Vector3d ray_normal = (ray_intersection-sphere_origin).normalized();

				// Simple diffuse model
				C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
		}
	}

	// Save to png
	write_matrix_to_png(C,C,C,A,filename);

}

void raytrace_sphere_perspective() {
	std::cout << "Simple ray tracer, one sphere with perspective projection" << std::endl;

	const std::string filename("sphere_perspective.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is perspecive, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d pix_origin(-1,1,0.9);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// Single light source
	const Vector3d light_position(-1,1,1);

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// Prepare the ray
			Vector3d ray_origin = origin;
			Vector3d ray_direction = (pix_origin + double(i)*x_displacement + double(j)*y_displacement - origin).normalized();
			
			// Intersect with the sphere
			Vector2d ray_on_xy(ray_origin(0),ray_origin(1));
			const Vector3d sphere_origin(0,0,0);
			const double sphere_radius = 0.9;
			Vector3d x = intersection_sphere(ray_origin, ray_direction, sphere_origin, sphere_radius);

			if (x(2) > 0) {
				// The ray hit the sphere, compute the exact intersection point
				Vector3d ray_intersection = ray_origin + std::min(x(0),x(1))*ray_direction;

				// Compute normal at the intersection point
				Vector3d ray_normal = (ray_intersection-sphere_origin).normalized();

				// Simple diffuse model
				C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
		}
	}

	// Save to png
	write_matrix_to_png(C,C,C,A,filename);

}

//function for the intersection
//e is the beginning of the vector, d is the direction
//a is the corner of the parallelogram, b-a and c-a are two sides
Vector3d intersection(Vector3d e, Vector3d d, Vector3d a, Vector3d b, Vector3d c) {
        MatrixXd T = MatrixXd::Zero(3,3);
        T << a-b, a-c, d;
        return T.colPivHouseholderQr().solve(a-e); 
    }

void raytrace_parallelogram() {
	std::cout << "Simple ray tracer, one parallelogram with orthographic projection" << std::endl;

	const std::string filename("plane_orthographic.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// TODO: Parameters of the parallelogram (position of the lower-left corner + two sides)
	Vector3d pgram_origin(0.2,-0.2,0.5);    
	Vector3d pgram_u(0,0.2,0.5);
	Vector3d pgram_v(-0.2,-0.2,0.5);

	// Single light source
	const Vector3d light_position(-1,1,1);

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// Prepare the ray
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = RowVector3d(0,0,-1);

			// TODO: Check if the ray intersects with the parallelogram
			Vector3d x = intersection(ray_origin, ray_direction, pgram_origin, pgram_u, pgram_v);
			if (x(0)>=0 && x(0)<=1 && x(1)>=0 && x(1)<=1 && x(2)>=0) {
				// TODO: The ray hit the parallelogram, compute the exact intersection point
				Vector3d ray_intersection = ray_origin + x(2)*ray_direction;

				// TODO: Compute normal at the intersection point
				Vector3d ray_normal = (pgram_u-pgram_origin).cross(pgram_v-pgram_origin).normalized();

				// Simple diffuse model
				C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
		}
	}

	// Save to png
	write_matrix_to_png(C,C,C,A,filename);
}

void raytrace_perspective() {
	std::cout << "Simple ray tracer, one parallelogram with perspective projection" << std::endl;

	const std::string filename("plane_perspective.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d pix_origin(-1,1,0.9);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// TODO: Parameters of the parallelogram (position of the lower-left corner + two sides)
	Vector3d pgram_origin(0.2,-0.2, 0.5);    
	Vector3d pgram_u(0,0.2, 0.5);
	Vector3d pgram_v(-0.2,-0.2, 0.5);

	// Single light source
	const Vector3d light_position(-1,1,1);

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// TODO: Prepare the ray (origin point and direction)
			Vector3d ray_origin = origin;
			Vector3d ray_direction = (pix_origin + double(i)*x_displacement + double(j)*y_displacement-origin).normalized();

			// TODO: Check if the ray intersects with the parallelogram
			Vector3d x = intersection(ray_origin, ray_direction, pgram_origin, pgram_u, pgram_v);
			if (x(0)>=0 && x(0)<=1 && x(1)>=0 && x(1)<=1 && x(2)>=0) {
				// TODO: The ray hit the parallelogram, compute the exact intersection point
				Vector3d ray_intersection = ray_origin + x(2)*ray_direction;

				// TODO: Compute normal at the intersection point
				Vector3d ray_normal = (pgram_u-pgram_origin).cross(pgram_v-pgram_origin).normalized();

				// Simple diffuse model
				C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
		}
	}

	// Save to png
	write_matrix_to_png(C,C,C,A,filename);
}

void raytrace_shading(){
	std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

	const std::string filename("shading.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d pix_origin(-1,1,0.9);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// Single light source
	const Vector3d light_position(-1,1,1);
	double ambient = 0.1;
	MatrixXd diffuse = MatrixXd::Zero(800, 800);
	MatrixXd specular = MatrixXd::Zero(800, 800);

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// Prepare the ray
			Vector3d ray_origin = origin;
			Vector3d ray_direction = (pix_origin + double(i)*x_displacement + double(j)*y_displacement - origin).normalized();
			
			// Intersect with the sphere
			Vector2d ray_on_xy(ray_origin(0),ray_origin(1));
			const Vector3d sphere_origin(0,0,0);
			const double sphere_radius = 0.9;
			Vector3d x = intersection_sphere(ray_origin, ray_direction, sphere_origin, sphere_radius);

			if (x(2) > 0) {
				// The ray hit the sphere, compute the exact intersection point
				Vector3d ray_intersection = ray_origin + std::min(x(0),x(1))*ray_direction;

				// Compute normal at the intersection point
				Vector3d ray_normal = (ray_intersection-sphere_origin).normalized();

				// TODO: Add shading parameter here
				double k_d = 2;
				double k_s = 0.5;

				diffuse(i,j) = k_d*(light_position-ray_intersection).normalized().transpose() * ray_normal;
				specular(i,j) = k_s*(origin - ray_intersection + light_position-ray_intersection).normalized().transpose() * ray_normal;

				// Simple diffuse model
				C(i,j) = ambient + diffuse(i,j) + specular(i,j);

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
		}
	}

	// Save to png
	write_matrix_to_png(C,C,C*0,A,filename);
}

int main() {
	raytrace_sphere();
	raytrace_sphere_perspective();
	raytrace_parallelogram();
	raytrace_perspective();
	raytrace_shading();
	
	
	return 0;
}
