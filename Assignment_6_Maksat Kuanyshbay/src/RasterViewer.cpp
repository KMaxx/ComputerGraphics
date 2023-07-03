#include "SDLViewer.h"


#include <Eigen/Core>

#include <functional>
#include <iostream>

#include "raster.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

Eigen::Vector2f intersection(Eigen::Vector4f p, Eigen::Vector4f a, Eigen::Vector4f b, Eigen::Vector4f c) {
    Eigen::MatrixXf T = Eigen::MatrixXf::Zero(2, 2);
    T << (a-b).head(2), (a-c).head(2);
    return T.inverse()*((a-p).head(2));
}

int factorial(int n) {
    if (n == 0) {
        return 1;
    }
    return n * factorial(n - 1);
}

int main(int argc, char *args[])
{
    int width = 500;
    int height = 500;
    // The Framebuffer storing the image rendered by the rasterizer
	Eigen::Matrix<FrameBufferAttributes,Eigen::Dynamic,Eigen::Dynamic> frameBuffer(width, height);
    
	// Global Constants (empty in this example)
	UniformAttributes uniform;
    uniform.clicks = 0;
    uniform.mode = 'd';
    uniform.is_pressed = false;
    uniform.is_animation = false;
    uniform.is_bezier = false;
    uniform.a = 0;
    uniform.b = 0;
    uniform.c = 0;
    uniform.col = 0;
    uniform.s << 0.0, 0.0;
    uniform.bc << 0.0,0.0,0.0,1.0;
    uniform.p << 0.0, 0.0, 0.0, 1.0;
    uniform.acolor << 1.0, 0.0, 0.0, 1.0;
    uniform.bcolor << 1.0, 0.0, 0.0, 1.0;
    uniform.ccolor << 1.0, 0.0, 0.0, 1.0;
    uniform.tc <<
        1.0, 0.0, 0.0, uniform.bc(0),
        0.0, 1.0, 0.0, uniform.bc(1),
        0.0, 0.0, 1.0, uniform.bc(2),
        0.0, 0.0, 0.0, 1.0;
    uniform.ct <<
        1.0, 0.0, 0.0, -uniform.bc(0),
        0.0, 1.0, 0.0, -uniform.bc(1),
        0.0, 0.0, 1.0, -uniform.bc(2),
        0.0, 0.0, 0.0, 1.0;
    uniform.rt <<
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0;
    uniform.view <<
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0;
    uniform.zoom <<
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0;
    Eigen::Matrix4f m_plus;
    m_plus <<
        1.2, 0.0, 0.0, 0.0,
        0.0, 1.2, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0;
    Eigen::Matrix4f m_minus;
    m_minus <<
        0.8, 0.0, 0.0, 0.0,
        0.0, 0.8, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0;
    Eigen::Matrix4f m_w;
    m_w <<
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, -0.4,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0;
    Eigen::Matrix4f m_a;
    m_a <<
        1.0, 0.0, 0.0, 0.4,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0;
    Eigen::Matrix4f m_s;
    m_s <<
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.4,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0;
    Eigen::Matrix4f m_d;
    m_d <<
        1.0, 0.0, 0.0, -0.4,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0;
    Eigen::Matrix4f m_t;
    m_t <<
        2.0/float(width), 0.0, 0.0, -1.0,
        0.0, -2.0/float(height), 0.0, 1.0-2.0/float(height),
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0;
	// Basic rasterization program
	Program program;

	// The vertex shader is the identity
	program.VertexShader = [](const VertexAttributes& va, const UniformAttributes& uniform)
	{
        VertexAttributes out;
        int n = uniform.trans.size();
        int in = va.index / 3;
        if (va.index > -1 && n > in) {
            out.position = uniform.zoom * uniform.trans[in] * va.position;
        }
        else {
            out.position = uniform.zoom * va.position;
        }
        
        out.color = va.color;
        out.index = va.index;
        return out;
	};

	// The fragment shader uses a fixed color
	program.FragmentShader = [](const VertexAttributes& va, const UniformAttributes& uniform)
	{
		return FragmentAttributes(va.color(0),va.color(1),va.color(2));
	};

	// The blending shader converts colors between 0 and 1 to uint8
	program.BlendingShader = [](const FragmentAttributes& fa, const FrameBufferAttributes& previous)
	{
		return FrameBufferAttributes(fa.color[0]*255,fa.color[1]*255,fa.color[2]*255,fa.color[3]*255);
	};

	// One triangle in the center of the screen
	std::vector<VertexAttributes> vertices;
    std::vector<VertexAttributes> lvertices;

    // Initialize the viewer and the corresponding callbacks
    SDLViewer viewer;
    viewer.init("Viewer Example", width, height);

    viewer.mouse_move = [&](int x, int y, int xrel, int yrel){
        Eigen::Vector4f temp(float(x), float(y), 0.0, 1.0);
        temp = m_t.inverse() * uniform.zoom.inverse() * m_t * temp;
        x = temp[0];
        y = temp[1];
        if (uniform.mode == 'i') {
            if (uniform.clicks == 1) {
                if (xrel != 0 || yrel != 0) {
                    lvertices.pop_back();
                    lvertices.push_back(VertexAttributes((float(x) / float(width) * 2.0) - 1.0, (float(height - 1.0 - y) / float(height) * 2.0) - 1.0, 0.0));
                    lvertices.back().color << 1.0, 0.0, 0.0, 1.0;
                    lvertices.back().index = -1;
                    viewer.redraw_next = true;
                }
            }
            else if (uniform.clicks == 2) {
                if (xrel != 0 || yrel != 0) {
                    vertices.pop_back();
                    vertices.push_back(VertexAttributes((float(x) / float(width) * 2.0) - 1.0, (float(height - 1.0 - y) / float(height) * 2.0) - 1.0, 0.0));
                    vertices.back().color << 1.0, 0.0, 0.0, 1.0;
                    vertices.back().index = vertices.size() - 1;
                    viewer.redraw_next = true;
                }
            }
        }
        else if (uniform.mode == 'o') {
            if (uniform.is_pressed) {
                if (xrel != 0 || yrel != 0) {
                    Eigen::Vector4f old_a = uniform.p;
                    uniform.p << (float(x) / float(width) * 2) - 1, (float(height - 1 - y) / float(height) * 2) - 1, 0, 1;
                    Eigen::Matrix4f m;
                    m <<
                        1.0, 0.0, 0.0, (uniform.p - old_a)[0],
                        0.0, 1.0, 0.0, (uniform.p - old_a)[1],
                        0.0, 0.0, 1.0, 0.0,
                        0.0, 0.0, 0.0, 1.0;
                    uniform.trans[uniform.a / 3] = m * uniform.trans[uniform.a / 3];
                    uniform.bc = m * uniform.bc;
                    uniform.tc <<
                        1.0, 0.0, 0.0, uniform.bc(0),
                        0.0, 1.0, 0.0, uniform.bc(1),
                        0.0, 0.0, 1.0, uniform.bc(2),
                        0.0, 0.0, 0.0, 1.0;
                    uniform.ct <<
                        1.0, 0.0, 0.0, -uniform.bc(0),
                        0.0, 1.0, 0.0, -uniform.bc(1),
                        0.0, 0.0, 1.0, -uniform.bc(2),
                        0.0, 0.0, 0.0, 1.0;
                    viewer.redraw_next = true;
                }
            }
        }
    };

    viewer.mouse_pressed = [&](int x, int y, bool is_pressed, int button, int clicks) {
        Eigen::Vector4f temp(float(x), float(y), 0.0, 1.0);
        temp = m_t.inverse() * uniform.zoom.inverse() * m_t * temp;
        x = temp[0];
        y = temp[1];
        if (uniform.mode == 'i') {
            if (is_pressed) {
                if (button == 1) {
                    if (clicks == 1) {
                        uniform.clicks = uniform.clicks + 1;
                        vertices.push_back(VertexAttributes((float(x) / float(width) * 2) - 1, (float(height - 1 - y) / float(height) * 2) - 1, 0));
                        vertices.back().color << 1, 0, 0, 1;
                        vertices.back().index = vertices.size()-1;
                        if (uniform.clicks == 1) {
                            lvertices.push_back(VertexAttributes((float(x) / float(width) * 2) - 1, (float(height - 1 - y) / float(height) * 2) - 1, 0));
                            lvertices.back().color << 1, 0, 0, 1;
                            lvertices.back().index = -1;
                            lvertices.push_back(VertexAttributes((float(x) / float(width) * 2) - 1, (float(height - 1 - y) / float(height) * 2) - 1, 0));
                            lvertices.back().color << 1, 0, 0, 1;
                            lvertices.back().index = -1;
                            viewer.redraw_next = true;
                        }
                        else if (uniform.clicks == 2) {
                            vertices.push_back(VertexAttributes((float(x) / float(width) * 2) - 1, (float(height - 1 - y) / float(height) * 2) - 1, 0));
                            vertices.back().color << 1, 0, 0, 1;
                            vertices.back().index = vertices.size() - 1;
                            lvertices.pop_back();
                            lvertices.pop_back();
                            viewer.redraw_next = true;
                        }
                        else if (uniform.clicks == 3) {
                            vertices.pop_back();
                            Eigen::Matrix4f m_e;
                            m_e <<
                                1.0, 0.0, 0.0, 0.0,
                                0.0, 1.0, 0.0, 0.0,
                                0.0, 0.0, 1.0, 0.0,
                                0.0, 0.0, 0.0, 1.0;
                            uniform.trans.push_back(m_e);
                            uniform.clicks = 0;
                            viewer.redraw_next = true;
                        }
                    }
                }
            }
        }
        else if (uniform.mode == 'o') {
            if (clicks == 1) {
                if (button == 1) {
                    if (is_pressed) {
                        uniform.p << (float(x) / float(width) * 2) - 1, (float(height - 1 - y) / float(height) * 2) - 1, 0, 1;
                        for (int i = 0;i < vertices.size()/3;i++) {
                            Eigen::Vector4f aa = uniform.trans[(vertices.size() - 3 - i * 3) / 3] * vertices[vertices.size() - 3 - i * 3].position;
                            Eigen::Vector4f bb = uniform.trans[(vertices.size() - 3 - i * 3) / 3] * vertices[vertices.size() - 2 - i * 3].position;
                            Eigen::Vector4f cc = uniform.trans[(vertices.size() - 3 - i * 3) / 3] * vertices[vertices.size() - 1 - i * 3].position;
                            uniform.s = intersection(uniform.p, aa, bb, cc);
                            if (uniform.s[0] > 0 && uniform.s[1] > 0 && uniform.s[0] + uniform.s[1] < 1) {
                                if (uniform.a != uniform.b){
                                    vertices[uniform.a].color = uniform.acolor;
                                    vertices[uniform.b].color = uniform.bcolor;
                                    vertices[uniform.c].color = uniform.ccolor;
                                }
                                uniform.a = vertices.size() - 3 - i * 3;
                                uniform.b = vertices.size() - 2 - i * 3;
                                uniform.c = vertices.size() - 1 - i * 3;
                                uniform.bc = (uniform.trans[uniform.a/3] * vertices[uniform.a].position + uniform.trans[uniform.a / 3] * vertices[uniform.b].position + uniform.trans[uniform.a / 3] * vertices[uniform.c].position) / 3;
                                uniform.tc <<
                                1.0, 0.0, 0.0, uniform.bc(0),
                                0.0, 1.0, 0.0, uniform.bc(1),
                                0.0, 0.0, 1.0, uniform.bc(2),
                                0.0, 0.0, 0.0, 1.0;
                                uniform.ct <<
                                1.0, 0.0, 0.0, -uniform.bc(0),
                                0.0, 1.0, 0.0, -uniform.bc(1),
                                0.0, 0.0, 1.0, -uniform.bc(2),
                                0.0, 0.0, 0.0, 1.0;
                                uniform.is_pressed = is_pressed;
                                uniform.acolor = vertices[uniform.a].color;
                                uniform.bcolor = vertices[uniform.b].color;
                                uniform.ccolor = vertices[uniform.c].color;
                                vertices[uniform.a].color << 0, 0, 1, 1;
                                vertices[uniform.b].color << 0, 0, 1, 1;
                                vertices[uniform.c].color << 0, 0, 1, 1;
                                viewer.redraw_next = true;
                                break; 
                            }
                        }
                    }
                    else {
                        uniform.is_pressed = false;
                        viewer.redraw_next = true;
                    }
                }
            }
        }
        else if (uniform.mode == 'p') {
            if (clicks == 1) {
                if (button == 1) {
                    if (is_pressed) {
                        VertexAttributes p = VertexAttributes((float(x) / float(width) * 2) - 1, (float(height - 1 - y) / float(height) * 2) - 1, 0);
                        p.color << 0, 0, 1, 1;
                        for (int i = 0;i < int(vertices.size()) / 3;i++) {
                            int n = vertices.size();
                            Eigen::Vector4f aa = uniform.trans[(n - 3 - i * 3) / 3] * vertices[n - 3 - i * 3].position;
                            Eigen::Vector4f bb = uniform.trans[(n - 3 - i * 3) / 3] * vertices[n - 2 - i * 3].position;
                            Eigen::Vector4f cc = uniform.trans[(n - 3 - i * 3) / 3] * vertices[n - 1 - i * 3].position;
                            uniform.s = intersection(p.position, aa, bb, cc);
                            if (uniform.s[0] > 0 && uniform.s[1] > 0 && uniform.s[0] + uniform.s[1] < 1) {
                                vertices.erase(vertices.begin() + n - 3 - i * 3, vertices.begin() + n - i * 3);
                                uniform.trans.erase(uniform.trans.begin() + (n - 3 - i * 3) / 3);
                                for (int j = n - 3 - i * 3; j < vertices.size(); j++) {
                                    vertices[j].index = j;
                                }
                                uniform.a = 0;
                                uniform.b = 0;
                                uniform.c = 0;
                                uniform.col = 0;
                                uniform.acolor << 1, 0, 0, 1;
                                uniform.bcolor << 1, 0, 0, 1;
                                uniform.ccolor << 1, 0, 0, 1;
                                uniform.s << 0, 0;
                                viewer.redraw_next = true;
                                break;
                            }
                        }
                    }
                }
            }

        }
        else if (uniform.mode == 'c') {
            if (clicks == 1) {
                if (button == 1) {
                    if (is_pressed) {
                        VertexAttributes p = VertexAttributes((float(x) / float(width) * 2) - 1, (float(height - 1 - y) / float(height) * 2) - 1, 0);
                        p.color << 0, 0, 1, 1;
                        float min = HUGE_VALF;
                        for (int i = 0;i < vertices.size();i++) {
                            if (min > (p.position - uniform.trans[(vertices.size() - 1 - i)/3]*vertices[vertices.size() - 1 - i].position).norm()) {
                                min = (p.position - uniform.trans[(vertices.size() - 1 - i) / 3] * vertices[vertices.size() - 1 - i].position).norm();
                                uniform.col = vertices.size() - 1 - i;
                            }
                        }
                    }
                }
            }
        }
    };

    viewer.mouse_wheel = [&](int dx, int dy, bool is_direction_normal) {
        for (int i = 0;i < vertices.size();i++) {
            std::cout << vertices[i].position;
        }
    };

    viewer.key_pressed = [&](char key, bool is_pressed, int modifier, int repeat) {
        if (is_pressed) {
            if (key == 'i' || key == 'o' || key == 'p' || key == 'c') {
                if (key != 'o' && uniform.mode == 'o' && vertices.size() != 0) {
                    vertices[uniform.a].color = uniform.acolor;
                    vertices[uniform.b].color = uniform.bcolor;
                    vertices[uniform.c].color = uniform.ccolor;
                    viewer.redraw_next = true;
                }
                if (uniform.mode == 'i' && vertices.size() != 0) {
                    if (vertices.size() % 3 == 1) {
                        lvertices.clear();
                        vertices.pop_back();
                    }
                    else if (vertices.size() % 3 == 2) {
                        lvertices.clear();
                        vertices.pop_back();
                        vertices.pop_back();
                    }
                }
                uniform.mode = key;
                uniform.clicks = 0;
            }
            else if (key == '=') {
                uniform.zoom = m_plus * uniform.zoom;
                viewer.redraw_next = true;
            }
            else if (key == '-') {
                uniform.zoom = m_minus * uniform.zoom;
                viewer.redraw_next = true;
            }
            else if (key == 'w') {
                uniform.zoom = m_w * uniform.zoom;
                viewer.redraw_next = true;
            }
            else if (key == 'a') {
                uniform.zoom = m_a * uniform.zoom;
                viewer.redraw_next = true;
            }
            else if (key == 's') {
                uniform.zoom = m_s * uniform.zoom;
                viewer.redraw_next = true;
            }
            else if (key == 'd') {
                uniform.zoom = m_d * uniform.zoom;
                viewer.redraw_next = true;
            }
            else if (uniform.mode == 'o') {
                if (key == 'm') {
                    uniform.keyfs.push_back(uniform.trans[uniform.a / 3]);
                }
                else if (key == 'x') {
                    uniform.keyfs.erase(uniform.keyfs.begin(), uniform.keyfs.end());
                }
                else if (key == 'n') {
                    uniform.is_animation = true;
                    viewer.redraw_next = true;
                }
                else if (key == 'b') {
                    uniform.is_bezier = true;
                    uniform.is_animation = true;
                    viewer.redraw_next = true;
                }
                else {
                    switch (key) {
                    case 'h':
                        uniform.rt <<
                            cos(0.174), sin(0.174), 0.0, 0.0,
                            -sin(0.174), cos(0.174), 0.0, 0.0,
                            0.0, 0.0, 1.0, 0.0,
                            0.0, 0.0, 0.0, 1.0;
                        break;
                    case 'j':
                        uniform.rt <<
                            cos(0.174), -sin(0.174), 0.0, 0.0,
                            sin(0.174), cos(0.174), 0.0, 0.0,
                            0.0, 0.0, 1.0, 0.0,
                            0.0, 0.0, 0.0, 1.0;
                        break;
                    case 'k':
                        uniform.rt <<
                            1.25, 0.0, 0.0, 0.0,
                            0.0, 1.25, 0.0, 0.0,
                            0.0, 0.0, 1.0, 0.0,
                            0.0, 0.0, 0.0, 1.0;
                        break;
                    case 'l':
                        uniform.rt <<
                            0.75, 0.0, 0.0, 0.0,
                            0.0, 0.75, 0.0, 0.0,
                            0.0, 0.0, 1.0, 0.0,
                            0.0, 0.0, 0.0, 1.0;
                        break;
                    default:
                        uniform.rt <<
                            1.0, 0.0, 0.0, 0.0,
                            0.0, 1.0, 0.0, 0.0,
                            0.0, 0.0, 1.0, 0.0,
                            0.0, 0.0, 0.0, 1.0;
                        break;
                    }
                    uniform.view = uniform.tc * uniform.rt * uniform.ct;
                    uniform.trans[uniform.a / 3] = uniform.view * uniform.trans[uniform.a / 3];
                    uniform.bc = uniform.view * uniform.bc;
                    uniform.tc <<
                        1.0, 0.0, 0.0, uniform.bc(0),
                        0.0, 1.0, 0.0, uniform.bc(1),
                        0.0, 0.0, 1.0, uniform.bc(2),
                        0.0, 0.0, 0.0, 1.0;
                    uniform.ct <<
                        1.0, 0.0, 0.0, -uniform.bc(0),
                        0.0, 1.0, 0.0, -uniform.bc(1),
                        0.0, 0.0, 1.0, -uniform.bc(2),
                        0.0, 0.0, 0.0, 1.0;
                    viewer.redraw_next = true;
                }
            }
            else if (uniform.mode == 'c') {
                switch (key) {
                case '1':
                    vertices[uniform.col].color << 0, 1, 0, 1;
                    break;
                case '2':
                    vertices[uniform.col].color << 0, 1, 1, 1;
                    break;
                case '3':
                    vertices[uniform.col].color << 1, 0, 1, 1;
                    break;
                case '4':
                    vertices[uniform.col].color << 0.5, 0.5, 0, 1;
                    break;
                case '5':
                    vertices[uniform.col].color << 0, 0.5, 0.5, 1;
                    break;
                case '6':
                    vertices[uniform.col].color << 1, 1, 1, 1;
                    break;
                case '7':
                    vertices[uniform.col].color << 0.5, 0.5, 0.5, 1;
                    break;
                case '8':
                    vertices[uniform.col].color << 0.5, 0, 0.5, 1;
                    break;
                case '9':
                    vertices[uniform.col].color << 0, 0, 0.5, 1;
                    break;
                default:
                    break;
                }
                if (key == '1' || key == '2' || key == '3' || key == '4' || key == '5' || key == '6' || key == '7' || key == '8' || key == '9') {
                    if (uniform.col == uniform.a) {
                        uniform.acolor = vertices[uniform.col].color;
                    }
                    if (uniform.col == uniform.b) {
                        uniform.bcolor = vertices[uniform.col].color;
                    }
                    if (uniform.col == uniform.c) {
                        uniform.ccolor = vertices[uniform.col].color;
                    }
                    viewer.redraw_next = true;
                }
            }
        }
    };

    viewer.redraw = [&](SDLViewer &viewer) {
        if (uniform.is_animation) {
            uniform.is_animation = false;
            if (uniform.is_bezier) {
                uniform.is_bezier = false;
                for (float t = 0.0; t < 1.1; t = t + 0.1) {
                    int n = uniform.keyfs.size();
                    uniform.trans[uniform.a / 3] = pow(1 - t, n-1) * uniform.keyfs[0];
                    for (int i = 1; i < n; i++) {
                        uniform.trans[uniform.a / 3] = uniform.trans[uniform.a / 3] + factorial(n-1) / (factorial(n-1 - i) * factorial(i)) * pow(1 - t, n-1 - i) * pow(t, i) * uniform.keyfs[i];
                    }
                    // Clear the framebuffer
                    for (unsigned i = 0;i < frameBuffer.rows();i++)
                        for (unsigned j = 0;j < frameBuffer.cols();j++)
                            frameBuffer(i, j).color << 0, 0, 0, 1;


                    rasterize_triangles(program, uniform, vertices, frameBuffer);
                    if (uniform.clicks == 1) {
                        rasterize_lines(program, uniform, lvertices, 1, frameBuffer);
                    }

                    // Buffer for exchanging data between rasterizer and sdl viewer
                    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> R(width, height);
                    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> G(width, height);
                    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> B(width, height);
                    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> A(width, height);

                    for (unsigned i = 0; i < frameBuffer.rows();i++)
                    {
                        for (unsigned j = 0; j < frameBuffer.cols();j++)
                        {
                            R(i, frameBuffer.cols() - 1 - j) = frameBuffer(i, j).color(0);
                            G(i, frameBuffer.cols() - 1 - j) = frameBuffer(i, j).color(1);
                            B(i, frameBuffer.cols() - 1 - j) = frameBuffer(i, j).color(2);
                            A(i, frameBuffer.cols() - 1 - j) = frameBuffer(i, j).color(3);
                        }
                    }
                    viewer.draw_image(R, G, B, A);
                }
            }
            else {
                for (int i = 0; i < uniform.keyfs.size() - 1; i++) {
                    for (float t = 0.0; t < 1.1; t = t + 0.1) {
                        uniform.trans[uniform.a / 3] = (1 - t) * uniform.keyfs[i] + t * uniform.keyfs[i + 1];
                        // Clear the framebuffer
                        for (unsigned i = 0;i < frameBuffer.rows();i++)
                            for (unsigned j = 0;j < frameBuffer.cols();j++)
                                frameBuffer(i, j).color << 0, 0, 0, 1;


                        rasterize_triangles(program, uniform, vertices, frameBuffer);
                        if (uniform.clicks == 1) {
                            rasterize_lines(program, uniform, lvertices, 1, frameBuffer);
                        }

                        // Buffer for exchanging data between rasterizer and sdl viewer
                        Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> R(width, height);
                        Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> G(width, height);
                        Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> B(width, height);
                        Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> A(width, height);

                        for (unsigned i = 0; i < frameBuffer.rows();i++)
                        {
                            for (unsigned j = 0; j < frameBuffer.cols();j++)
                            {
                                R(i, frameBuffer.cols() - 1 - j) = frameBuffer(i, j).color(0);
                                G(i, frameBuffer.cols() - 1 - j) = frameBuffer(i, j).color(1);
                                B(i, frameBuffer.cols() - 1 - j) = frameBuffer(i, j).color(2);
                                A(i, frameBuffer.cols() - 1 - j) = frameBuffer(i, j).color(3);
                            }
                        }
                        viewer.draw_image(R, G, B, A);
                    }
                }
            }
        }
        else {
            // Clear the framebuffer
            for (unsigned i = 0;i < frameBuffer.rows();i++)
                for (unsigned j = 0;j < frameBuffer.cols();j++)
                    frameBuffer(i, j).color << 0, 0, 0, 1;


            rasterize_triangles(program, uniform, vertices, frameBuffer);
            if (uniform.clicks == 1) {
                rasterize_lines(program, uniform, lvertices, 1, frameBuffer);
            }

            // Buffer for exchanging data between rasterizer and sdl viewer
            Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> R(width, height);
            Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> G(width, height);
            Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> B(width, height);
            Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> A(width, height);

            for (unsigned i = 0; i < frameBuffer.rows();i++)
            {
                for (unsigned j = 0; j < frameBuffer.cols();j++)
                {
                    R(i, frameBuffer.cols() - 1 - j) = frameBuffer(i, j).color(0);
                    G(i, frameBuffer.cols() - 1 - j) = frameBuffer(i, j).color(1);
                    B(i, frameBuffer.cols() - 1 - j) = frameBuffer(i, j).color(2);
                    A(i, frameBuffer.cols() - 1 - j) = frameBuffer(i, j).color(3);
                }
            }
            viewer.draw_image(R, G, B, A);
        }
        
    };

    viewer.launch();

    return 0;
}
