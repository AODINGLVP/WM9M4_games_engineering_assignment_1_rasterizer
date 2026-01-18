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
        int W = renderer.canvas.getWidth();
        int H = renderer.canvas.getHeight();

        vec2D minV, maxV;
        getBoundsWindow(renderer.canvas, minV, maxV);

        int minX = clamp((int)std::floor(minV.x), 0, W - 1);
        int minY = clamp((int)std::floor(minV.y), 0, H - 1);
        int maxX = clamp((int)std::ceil(maxV.x), 0, W - 1);
        int maxY = clamp((int)std::ceil(maxV.y), 0, H - 1);

                   
        if (area < 1.f) return;
       
        float invArea = 1.0f / area;

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
        float w2_stepx = e2.x, w2_stepy = e2.y;
        float alpha;
        float beta ;
        float gamma;
        L.omega_i.normalise(); // 只做一次

        for (int y = minY; y <= maxY; ++y) {
            float w0 = row_w0, w1 = row_w1, w2 = row_w2;

            for (int x = minX; x <= maxX; ++x) {
                if (w0 >= 0.f && w1 >= 0.f && w2 >= 0.f) {
                     alpha = w0 * invArea;
                     beta = w1 * invArea;
                     gamma = w2 * invArea;

                    float depth = interpolate(alpha, beta, gamma, v[0].p[2], v[1].p[2], v[2].p[2]);

                    if (renderer.zbuffer(x, y) > depth && depth > 0.001f) {
                        colour c = interpolate(alpha, beta, gamma, v[0].rgb, v[1].rgb, v[2].rgb);
                        c.clampColour();
                        vec4 normal = interpolate(alpha, beta, gamma, v[0].normal, v[1].normal, v[2].normal);
                        normal.normalise();

                        float dot = std::max(vec4::dot(L.omega_i, normal), 0.0f);
                        colour a = (c * kd) * (L.L * dot) + (L.ambient * ka);

                        unsigned char r, g, b;
                        a.toRGB(r, g, b);
                        renderer.canvas.draw(x, y, r, g, b);
                        renderer.zbuffer(x, y) = depth;
                    }
                }
                w0 += w0_stepx; w1 += w1_stepx; w2 += w2_stepx;
            }
            row_w0 += w0_stepy; row_w1 += w1_stepy; row_w2 += w2_stepy;
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
