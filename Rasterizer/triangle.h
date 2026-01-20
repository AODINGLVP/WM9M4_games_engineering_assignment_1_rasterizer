#pragma once

#include "mesh.h"
#include "colour.h"
#include "renderer.h"
#include "light.h"
#include <iostream>
#include <algorithm>
#include <cmath>

// Simple support class for a 2D vector
class vec2D {
public:
    float x, y;

    // Default constructor initializes both components to 0
    vec2D() { x = y = 0.f; };

    // Constructor initializes components with given values
    vec2D(float _x, float _y) : x(_x), y(_y) {}

    // Constructor initializes components from a vec4
    vec2D(vec4 v) {
        x = v[0];
        y = v[1];
    }

    // Display the vector components
    void display() { std::cout << x << '\t' << y << std::endl; }

    // Overloaded subtraction operator for vector subtraction
    vec2D operator- (vec2D& v) {
        vec2D q;
        q.x = x - v.x;
        q.y = y - v.y;
        return q;
    }
};

// Class representing a triangle for rendering purposes
class triangle {
    Vertex v[3];       // Vertices of the triangle
    float area;        // Area of the triangle
    colour col[3];     // Colors for each vertex of the triangle

public:
    // Constructor initializes the triangle with three vertices
    // Input Variables:
    // - v1, v2, v3: Vertices defining the triangle
    triangle(const Vertex& v1, const Vertex& v2, const Vertex& v3) {
        v[0] = v1;
        v[1] = v2;
        v[2] = v3;

        // Calculate the 2D area of the triangle
        vec2D e1 = vec2D(v[1].p - v[0].p);
        vec2D e2 = vec2D(v[2].p - v[0].p);
        area = std::fabs(e1.x * e2.y - e1.y * e2.x);
    }

    // Helper function to compute the cross product for barycentric coordinates
    // Input Variables:
    // - v1, v2: Edges defining the vector
    // - p: Point for which coordinates are being calculated
    float getC(vec2D v1, vec2D v2, vec2D p) {
        vec2D e = v2 - v1;
        vec2D q = p - v1;
        return q.y * e.x - q.x * e.y;
    }

    // Compute barycentric coordinates for a given point
    // Input Variables:
    // - p: Point to check within the triangle
    // Output Variables:
    // - alpha, beta, gamma: Barycentric coordinates of the point
    // Returns true if the point is inside the triangle, false otherwise
    bool getCoordinates(vec2D p, float& alpha, float& beta, float& gamma) {
        alpha = getC(vec2D(v[0].p), vec2D(v[1].p), p) / area;
        beta = getC(vec2D(v[1].p), vec2D(v[2].p), p) / area;
        gamma = getC(vec2D(v[2].p), vec2D(v[0].p), p) / area;

        if (alpha < 0.f || beta < 0.f || gamma < 0.f) return false;
        return true;
    }

    // Template function to interpolate values using barycentric coordinates
    // Input Variables:
    // - alpha, beta, gamma: Barycentric coordinates
    // - a1, a2, a3: Values to interpolate
    // Returns the interpolated value
    template <typename T>
    T interpolate(float alpha, float beta, float gamma, T a1, T a2, T a3) {
        return (a1 * alpha) + (a2 * beta) + (a3 * gamma);
    }
    template <typename T>
    T compute_incremental_interpolate(float x0, float x1, float y0, float y1, float x2, float y2, T v0, T v1, T v2, float invarea, T& dvdx, T& dvdy) {
        dvdx = (v0 * (y1 - y2) + v1 * (y2 - y0) + v2 * (y0 - y1)) * invarea;
        dvdy = (v0 * (x2 - x1) + v1 * (x0 - x2) + v2 * (x1 - x0)) * invarea;
        return dvdx;
    }
    static vec4 makeEdge(const vec4& a, const vec4& b) {
        vec4 e(0,0,0,0);
        e.x = a.y - b.y;
        e.y = b.x - a.x;
        e.z = a.x * b.y - a.y * b.x;
        return e;
    }
    static  float evalEdge(const vec4& e, float x, float y) {
        return e.x * x + e.y * y + e.z;
    }
    // Draw the triangle on the canvas
    // Input Variables:
    // - renderer: Renderer object for drawing
    // - L: Light object for shading calculations
    // - ka, kd: Ambient and diffuse lighting coefficients
    void draw(Renderer& renderer, Light& L, float ka, float kd) {


        if (area < 1.f) return;
       


        float invArea = 1.0f / area;
      

        int W = renderer.canvas.getWidth();
        int H = renderer.canvas.getHeight();

        vec2D minV, maxV;
        getBoundsWindow(renderer.canvas, minV, maxV);

        int minX = clamp((int)std::floor(minV.x), 0, W - 1);
        int minY = clamp((int)std::floor(minV.y), 0, H - 1);
        int maxX = clamp((int)std::ceil(maxV.x), 0, W - 1);
        int maxY = clamp((int)std::ceil(maxV.y), 0, H - 1);

    
		float dzdx, dzdy;
		compute_incremental_interpolate(v[0].p.x, v[1].p.x, v[0].p.y, v[1].p.y, v[2].p.x, v[2].p.y, v[0].p.z, v[1].p.z, v[2].p.z, invArea, dzdx, dzdy);
		colour dcdx, dcdy;
		compute_incremental_interpolate(v[0].p.x, v[1].p.x, v[0].p.y, v[1].p.y, v[2].p.x, v[2].p.y, v[0].rgb, v[1].rgb, v[2].rgb, invArea, dcdx, dcdy);
		vec4 dndx, dndy;
		compute_incremental_interpolate(v[0].p.x, v[1].p.x, v[0].p.y, v[1].p.y, v[2].p.x, v[2].p.y, v[0].normal, v[1].normal, v[2].normal, invArea, dndx, dndy);

        vec4 e0 = makeEdge(v[1].p, v[2].p); // (A,B,C)
        vec4 e1 = makeEdge(v[2].p, v[0].p);
        vec4 e2 = makeEdge(v[0].p, v[1].p);

        float startX = (float)minX;
        float startY = (float)minY;

        float row_w0 = evalEdge(e0, startX, startY);
        float row_w1 = evalEdge(e1, startX, startY);
        float row_w2 = evalEdge(e2, startX, startY);

        float w0_stepx = e0.x, w0_stepy = e0.y;
        float w1_stepx = e1.x, w1_stepy = e1.y;
        float w2_stepx = e2.x , w2_stepy = e2.y;
        float w0_block = e0.x * 8;      
        float w1_block = e1.x * 8;
        float w2_block = e2.x * 8;

		float alpha, beta, gamma;
        float alpha0 = row_w0 * invArea;
        float beta0 = row_w1 * invArea;
        float gamma0 = row_w2 * invArea;
		//std::cout << "asdasdds"<< alpha0+ beta0 + gamma0 <<std::endl;
        float z_row = v[0].p.z * alpha0 + v[1].p.z * beta0 + v[2].p.z * gamma0;
        colour c_row = v[0].rgb * alpha0 + v[1].rgb * beta0 + v[2].rgb * gamma0;
        vec4 n_row = v[0].normal * alpha0 + v[1].normal * beta0 + v[2].normal * gamma0;



        L.omega_i.normalise(); // 只做一次
        __m256 zero = _mm256_setzero_ps();
        alignas(32)float w0[8] = { row_w0 ,row_w0 + w0_stepx * 1,row_w0 + w0_stepx * 2,row_w0 + w0_stepx * 3,row_w0 + w0_stepx * 4,row_w0 + w0_stepx * 5,row_w0 + w0_stepx * 6,row_w0 + w0_stepx * 7 };
        alignas(32)float w1[8] = { row_w1 ,row_w1 + w1_stepx * 1,row_w1 + w1_stepx * 2,row_w1 + w1_stepx * 3,row_w1 + w1_stepx * 4,row_w1 + w1_stepx * 5,row_w1 + w1_stepx * 6,row_w1 + w1_stepx * 7 };
        alignas(32) float w2[8] = { row_w2 ,row_w2 + w2_stepx * 1,row_w2 + w2_stepx * 2,row_w2 + w2_stepx * 3,row_w2 + w2_stepx * 4,row_w2 + w2_stepx * 5,row_w2 + w2_stepx * 6,row_w2 + w2_stepx * 7 };
        


        for (int y = minY; y <= maxY; ++y) {
            for (int i = 0; i < 8; i++) {
                w0[i] = row_w0 + w0_stepx * i;
                w1[i] = row_w1 + w1_stepx * i;
				w2[i] = row_w2 + w2_stepx * i;
            }

            __m256 w0_vec = _mm256_load_ps(w0);
            __m256 w1_vec = _mm256_load_ps(w1);
            __m256 w2_vec = _mm256_load_ps(w2);
            __m256 w0_stepx_v = _mm256_set1_ps(w0_block);
            __m256 w1_stepx_v = _mm256_set1_ps(w1_block);
            __m256 w2_stepx_v = _mm256_set1_ps(w2_block);


            colour col = c_row;
            float z = z_row;
            vec4 nor = n_row;
           
            for (int x = minX; x <= maxX; x+=8) {
               
               
                __m256 m0 = _mm256_cmp_ps(w0_vec, zero, _CMP_GE_OQ);
                __m256 m1 = _mm256_cmp_ps(w1_vec, zero, _CMP_GE_OQ);
                __m256 m2 = _mm256_cmp_ps(w2_vec, zero, _CMP_GE_OQ);
                __m256 inside = _mm256_and_ps(_mm256_and_ps(m0, m1), m2);
                int mask = _mm256_movemask_ps(inside);
                if (mask == 0) {
                    w0_vec = _mm256_add_ps(w0_vec, w0_stepx_v);
                    w1_vec = _mm256_add_ps(w1_vec, w1_stepx_v);
                    w2_vec = _mm256_add_ps(w2_vec, w2_stepx_v);
                    z += dzdx * 8;
                    col = col + dcdx * 8;
                    nor = nor + dndx * 8;
                    continue;
                }

               
                while (mask) {
                    int i = _tzcnt_u32(mask);
                    mask &= (mask - 1);
                    

                   
                   

                     
                        float depth = z+ dzdx * i;
                      

                        if (renderer.zbuffer(x+i, y) > depth && depth > 0.001f) {
                         
                            colour c = col + dcdx * i;
                            c.clampColour();
							vec4 normal = nor+ dndx * i;
                            normal.normalise();
                            float dot = std::max(vec4::dot(L.omega_i, normal), 0.0f);
                            colour a = (c * kd) * (L.L * dot) + (L.ambient * ka);
                            unsigned char r, g, b;
                            a.toRGB(r, g, b);
                            renderer.canvas.draw(x+i, y, r, g, b);
                            renderer.zbuffer(x+i, y) = depth;
                        }
                    
                }
					w0_vec = _mm256_add_ps(w0_vec, w0_stepx_v);
					w1_vec = _mm256_add_ps(w1_vec, w1_stepx_v);
					w2_vec = _mm256_add_ps(w2_vec, w2_stepx_v);
                    z += dzdx*8;
                    col = col + dcdx * 8;
                    nor = nor + dndx * 8;
                    
                    if (maxX - x < 8) {
                        _mm256_store_ps(w0, w0_vec);
                        _mm256_store_ps(w1, w1_vec);
                        _mm256_store_ps(w2, w2_vec);
                        for (int i = 0; i <= maxX-x; i++) {
                            if (w0[i] >=0 && w1[i] >= 0 && w2[i] >= 0) {
                                float depth = z + dzdx * i;
                            

                                if (renderer.zbuffer(x + i, y) > depth && depth > 0.001f) {
                                    
                                    colour c = col + dcdx * i;
                                    c.clampColour();
                                    vec4 normal = nor + dndx * i;
                                    normal.normalise();
                                    float dot = std::max(vec4::dot(L.omega_i, normal), 0.0f);
                                    colour a = (c * kd) * (L.L * dot) + (L.ambient * ka);
                                    unsigned char r, g, b;
                                    a.toRGB(r, g, b);
                                    renderer.canvas.draw(x + i, y, r, g, b);
                                    renderer.zbuffer(x + i, y) = depth;
                                }
                            }
                        }
                        break;

                }
                
            }
            row_w0 += w0_stepy; row_w1 += w1_stepy; row_w2 += w2_stepy;
            z_row += dzdy;
            c_row = c_row + dcdy;
            n_row = n_row + dndy;
        }
    }


    // Compute the 2D bounds of the triangle
    // Output Variables:
    // - minV, maxV: Minimum and maximum bounds in 2D space
    void getBounds(vec2D& minV, vec2D& maxV) {
        minV = vec2D(v[0].p);
        maxV = vec2D(v[0].p);
        for (unsigned int i = 1; i < 3; i++) {
            minV.x = std::min(minV.x, v[i].p[0]);
            minV.y = std::min(minV.y, v[i].p[1]);
            maxV.x = std::max(maxV.x, v[i].p[0]);
            maxV.y = std::max(maxV.y, v[i].p[1]);
        }
    }

    // Compute the 2D bounds of the triangle, clipped to the canvas
    // Input Variables:
    // - canvas: Reference to the rendering canvas
    // Output Variables:
    // - minV, maxV: Clipped minimum and maximum bounds
    void getBoundsWindow(GamesEngineeringBase::Window& canvas, vec2D& minV, vec2D& maxV) {
        getBounds(minV, maxV);
        minV.x = std::max(minV.x, static_cast<float>(0));
        minV.y = std::max(minV.y, static_cast<float>(0));
        maxV.x = std::min(maxV.x, static_cast<float>(canvas.getWidth()));
        maxV.y = std::min(maxV.y, static_cast<float>(canvas.getHeight()));
    }

    // Debugging utility to display the triangle bounds on the canvas
    // Input Variables:
    // - canvas: Reference to the rendering canvas
    void drawBounds(GamesEngineeringBase::Window& canvas) {
        vec2D minV, maxV;
        getBounds(minV, maxV);

        for (int y = (int)minV.y; y < (int)maxV.y; y++) {
            for (int x = (int)minV.x; x < (int)maxV.x; x++) {
                canvas.draw(x, y, 255, 0, 0);
            }
        }
    }

    // Debugging utility to display the coordinates of the triangle vertices
    void display() {
        for (unsigned int i = 0; i < 3; i++) {
            v[i].p.display();
        }
        std::cout << std::endl;
    }
};
