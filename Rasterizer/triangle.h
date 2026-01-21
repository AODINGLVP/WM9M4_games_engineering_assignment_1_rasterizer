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


        alignas(32)float w0_step_vx[8] = { 0 ,  w0_stepx * 1, + w0_stepx * 2, + w0_stepx * 3, + w0_stepx * 4, + w0_stepx * 5, + w0_stepx * 6, + w0_stepx * 7 };
        alignas(32)float w1_step_vx[8] = { 0 ,  w1_stepx * 1, + w1_stepx * 2, + w1_stepx * 3, + w1_stepx * 4, + w1_stepx * 5, + w1_stepx * 6, + w1_stepx * 7 };
        alignas(32) float w2_step_vx[8] = { 0 ,  w2_stepx * 1, + w2_stepx * 2, + w2_stepx * 3, + w2_stepx * 4, + w2_stepx * 5, + w2_stepx * 6, + w2_stepx * 7 };
        __m256 w0_step_vx_256 = _mm256_load_ps(w0_step_vx);
        __m256 w1_step_vx_256 = _mm256_load_ps(w1_step_vx);
        __m256 w2_step_vx_256 = _mm256_load_ps(w2_step_vx);





        int zrow;
        __m256 w0_stepx_v = _mm256_set1_ps(w0_block);
        __m256 w1_stepx_v = _mm256_set1_ps(w1_block);
        __m256 w2_stepx_v = _mm256_set1_ps(w2_block);

        // for color
		//__m256 col_v_r = _mm256_set1_ps(c_row.r);
		//__m256 col_v_g = _mm256_set1_ps(c_row.g);
		//__m256 col_v_b = _mm256_set1_ps(c_row.b);
        

        alignas(32)float dcdx_step_vr[8] = { 0 ,  dndx.x * 1, +dndx.x * 2, +dndx.x * 3, +dndx.x * 4, +dndx.x * 5, +dndx.x * 6, +dndx.x * 7 };
        alignas(32)float dcdx_step_vg[8] = { 0 ,  dndx.y * 1, +dndx.y * 2, +dndx.y * 3, +dndx.y * 4, +dndx.y * 5, +dndx.y * 6, +dndx.y * 7 };
        alignas(32) float dcdx_step_vb[8] = { 0 ,  dndx.z * 1, +dndx.z * 2, +dndx.z * 3, +dndx.z * 4, +dndx.z * 5, +dndx.z * 6, +dndx.z * 7 };
		__m256 dcdx_v_r = _mm256_load_ps(dcdx_step_vr);
		__m256 dcdx_v_g = _mm256_load_ps(dcdx_step_vg);
        __m256 dcdx_v_b = _mm256_load_ps(dcdx_step_vb);


		colour dcdx_block = dcdx * 8;

        __m256 dcdx_block_vr = _mm256_set1_ps(dcdx_block.r);
        __m256 dcdx_block_vg = _mm256_set1_ps(dcdx_block.g);
        __m256 dcdx_block_vb = _mm256_set1_ps(dcdx_block.b);



        // for depth
       // __m256 depth_v = _mm256_set1_ps(z_row);
        
        alignas(32)float dzdx_step_vr[8] = { 0 ,  dzdx * 1, +dzdx * 2, +dzdx * 3, +dzdx * 4, +dzdx * 5, +dzdx * 6, +dzdx * 7 };
        __m256 dzdx_v_r = _mm256_load_ps(dzdx_step_vr);
       
        __m256 dzdx_block_v = _mm256_set1_ps(dzdx*8);
       
        //for normal
        alignas(32)float dndx_step_vx[8] = { 0 ,  dndx.x * 1, +dndx.x * 2, +dndx.x * 3, +dndx.x * 4, +dndx.x * 5, +dndx.x * 6, +dndx.x * 7 };
        alignas(32)float dndx_step_vy[8] = { 0 ,  dndx.y * 1, +dndx.y * 2, +dndx.y * 3, +dndx.y * 4, +dndx.y * 5, +dndx.y * 6, +dndx.y * 7 };
        alignas(32) float dndx_step_vz[8] = { 0 ,  dndx.z * 1, +dndx.z * 2, +dndx.z * 3, +dndx.z * 4, +dndx.z * 5, +dndx.z * 6, +dndx.z * 7 };
        __m256 dndx_v_x = _mm256_load_ps(dndx_step_vx);
        __m256 dndx_v_y = _mm256_load_ps(dndx_step_vy);
        __m256 dndx_v_z = _mm256_load_ps(dndx_step_vz);


        vec4 dndx_block = dndx * 8;

        __m256 dndx_block_vx = _mm256_set1_ps(dndx_block.x);
        __m256 dndx_block_vy = _mm256_set1_ps(dndx_block.y);
        __m256 dndx_block_vz = _mm256_set1_ps(dndx_block.z);



        //for light
        __m256 light_omega_vx = _mm256_set1_ps(L.omega_i.x);
        __m256 light_omega_vy = _mm256_set1_ps(L.omega_i.y);
        __m256 light_omega_vz = _mm256_set1_ps(L.omega_i.z);

        __m256 light_l_vr = _mm256_set1_ps(L.L.r);
        __m256 light_l_vg = _mm256_set1_ps(L.L.g);
        __m256 light_l_vb = _mm256_set1_ps(L.L.b);

        __m256 light_ambient_vr = _mm256_set1_ps(L.ambient.r);
        __m256 light_ambient_vg = _mm256_set1_ps(L.ambient.g);
        __m256 light_ambient_vb = _mm256_set1_ps(L.ambient.b);

        __m256 light_ka = _mm256_set1_ps(ka);
        __m256 light_kd = _mm256_set1_ps(kd);

		__m256 one = _mm256_set1_ps(1.0f);
        __m256 zeor_dot_oneoneone = _mm256_set1_ps(0.001f);
		__m256 twofivefive = _mm256_set1_ps(255.0f);

        for (int y = minY; y <= maxY; ++y) {
           
            __m256 w0_vec = _mm256_set1_ps(row_w0);
            __m256 w1_vec = _mm256_set1_ps(row_w1);
            __m256 w2_vec = _mm256_set1_ps(row_w2);
            w0_vec = _mm256_add_ps(w0_vec, w0_step_vx_256);
             w1_vec = _mm256_add_ps(w1_vec, w1_step_vx_256);
             w2_vec = _mm256_add_ps(w2_vec, w2_step_vx_256);
            
             //color
             __m256 col_v_r = _mm256_set1_ps(c_row.r);
             __m256 col_v_g = _mm256_set1_ps(c_row.g);
             __m256 col_v_b = _mm256_set1_ps(c_row.b);

             col_v_r = _mm256_add_ps(col_v_r, dcdx_v_r);
              col_v_g = _mm256_add_ps(col_v_g, dcdx_v_g);
              col_v_b = _mm256_add_ps(col_v_b, dcdx_v_b);
              //depth
			  __m256 depth_v = _mm256_set1_ps(z_row);
			  depth_v = _mm256_add_ps(depth_v, dzdx_v_r);

              //normal
              __m256 normal_v_x = _mm256_set1_ps(n_row.x);
              __m256 normal_v_y = _mm256_set1_ps(n_row.y);
              __m256 normal_v_z = _mm256_set1_ps(n_row.z);

              normal_v_x = _mm256_add_ps(normal_v_x, dndx_v_x);
              normal_v_y = _mm256_add_ps(normal_v_y, dndx_v_y);
              normal_v_z = _mm256_add_ps(normal_v_z, dndx_v_z);


            //colour col = c_row;

           // float z = z_row;
            //vec4 nor = n_row;
           
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
                    //++dcdx
                    col_v_r = _mm256_add_ps(col_v_r, dcdx_block_vr);
                    col_v_g = _mm256_add_ps(col_v_g, dcdx_block_vg);
                    col_v_b = _mm256_add_ps(col_v_b, dcdx_block_vb);
                    //++dzdx
					depth_v = _mm256_add_ps(depth_v, dzdx_block_v);
					//++dndx
                    normal_v_x = _mm256_add_ps(normal_v_x, dndx_block_vx);
                    normal_v_y = _mm256_add_ps(normal_v_y, dndx_block_vy);
                    normal_v_z = _mm256_add_ps(normal_v_z, dndx_block_vz);


                   // z += dzdx * 8;
                    //col = col + dcdx * 8;
                   // nor = nor + dndx * 8;
                    continue;
                }
                /*
                {
                    //colour c = col + dcdx * i;

                    //c.clampColour();

                    col_v_r = _mm256_min_ps(col_v_r, one);
                    //vec4 normal = nor + dndx * i;
                   
                   // normal.normalise();
                    __m256 normal_length = _mm256_sqrt_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(normal_v_x, normal_v_x),_mm256_mul_ps(normal_v_y, normal_v_y)),_mm256_mul_ps(normal_v_z, normal_v_z)));
                    normal_v_x = _mm256_div_ps(normal_v_x, normal_length);
                    normal_v_y = _mm256_div_ps(normal_v_y, normal_length);
                    normal_v_z = _mm256_div_ps(normal_v_z, normal_length); 
                   // float dot = std::max(vec4::dot(L.omega_i, normal), 0.0f);
					__m256 dot_v = _mm256_max_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(normal_v_x, light_omega_vx), _mm256_mul_ps(normal_v_y, light_omega_vy)), _mm256_mul_ps(normal_v_z, light_omega_vz)), zero);
                    
                   // colour a = (c * kd) * (L.L * dot) + (L.ambient * ka);
                    __m256 a_vr = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(col_v_r, light_kd), _mm256_mul_ps(light_l_vr, dot_v)), _mm256_mul_ps(light_ambient_vr, light_ka));
                    __m256 a_vg = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(col_v_g, light_kd), _mm256_mul_ps(light_l_vg, dot_v)), _mm256_mul_ps(light_ambient_vg, light_ka));
                    __m256 a_vb = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(col_v_b, light_kd), _mm256_mul_ps(light_l_vb, dot_v)), _mm256_mul_ps(light_ambient_vb, light_ka));


                    alignas(32) float a_r[8];
                    alignas(32) float a_g[8];
                    alignas(32) float a_b[8];
                    _mm256_store_ps(a_r, a_vr);
                    _mm256_store_ps(a_g, a_vg);
                    _mm256_store_ps(a_b, a_vb);
                    //extract the depth data from the SIMD register
                    alignas(32) float depth_8f[8];
                    _mm256_store_ps(depth_8f, depth_v);
                    unsigned char r, g, b;
                    float depth_end;
                    //可以考虑顺序使用八个r,八个g,八个b
                    for(int i=0;i<8;i++)
                    {
						r = a_r[i] * 255;
						g = a_g[i] * 255;
                        b = a_b[i] * 255;
                        depth_end = depth_8f[i];
                       
                            renderer.canvas.draw(x + i, y, r, g, b);
                            renderer.zbuffer(x + i, y) = depth_end;
                          
						
					}
                    

                    
                   
                   
                }
               */
				__m256 zbuffer_v = _mm256_loadu_ps(&renderer.zbuffer(x, y));
                __m256 m0_depth = _mm256_cmp_ps(depth_v, zeor_dot_oneoneone, _CMP_GE_OQ);
                __m256 m1_depth = _mm256_cmp_ps(zbuffer_v, depth_v, _CMP_GE_OQ);
                __m256 inside_depth = _mm256_and_ps(m0_depth, m1_depth);
                int mask_depth = _mm256_movemask_ps(inside_depth);

                if (mask_depth == 0) {
                    w0_vec = _mm256_add_ps(w0_vec, w0_stepx_v);
                    w1_vec = _mm256_add_ps(w1_vec, w1_stepx_v);
                    w2_vec = _mm256_add_ps(w2_vec, w2_stepx_v);
                    //++dcdx
                    col_v_r = _mm256_add_ps(col_v_r, dcdx_block_vr);
                    col_v_g = _mm256_add_ps(col_v_g, dcdx_block_vg);
                    col_v_b = _mm256_add_ps(col_v_b, dcdx_block_vb);
                    //++dzdx
                    depth_v = _mm256_add_ps(depth_v, dzdx_block_v);
                    //++dndx
                    normal_v_x = _mm256_add_ps(normal_v_x, dndx_block_vx);
                    normal_v_y = _mm256_add_ps(normal_v_y, dndx_block_vy);
                    normal_v_z = _mm256_add_ps(normal_v_z, dndx_block_vz);


                    // z += dzdx * 8;
                     //col = col + dcdx * 8;
                    // nor = nor + dndx * 8;
                    continue;
                }
				
				
               
                  
                         
                            {
                                //colour c = col + dcdx * i;

                                //c.clampColour();

                                col_v_r = _mm256_min_ps(col_v_r, one);
                                //vec4 normal = nor + dndx * i;

                               // normal.normalise();
                                __m256 normal_length = _mm256_sqrt_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(normal_v_x, normal_v_x), _mm256_mul_ps(normal_v_y, normal_v_y)), _mm256_mul_ps(normal_v_z, normal_v_z)));
                                normal_v_x = _mm256_div_ps(normal_v_x, normal_length);
                                normal_v_y = _mm256_div_ps(normal_v_y, normal_length);
                                normal_v_z = _mm256_div_ps(normal_v_z, normal_length);
                                // float dot = std::max(vec4::dot(L.omega_i, normal), 0.0f);
                                __m256 dot_v = _mm256_max_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(normal_v_x, light_omega_vx), _mm256_mul_ps(normal_v_y, light_omega_vy)), _mm256_mul_ps(normal_v_z, light_omega_vz)), zero);

                                // colour a = (c * kd) * (L.L * dot) + (L.ambient * ka);
                                __m256 a_vr = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(col_v_r, light_kd), _mm256_mul_ps(light_l_vr, dot_v)), _mm256_mul_ps(light_ambient_vr, light_ka));
                                __m256 a_vg = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(col_v_g, light_kd), _mm256_mul_ps(light_l_vg, dot_v)), _mm256_mul_ps(light_ambient_vg, light_ka));
                                __m256 a_vb = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(col_v_b, light_kd), _mm256_mul_ps(light_l_vb, dot_v)), _mm256_mul_ps(light_ambient_vb, light_ka));


                                int realmask = _mm256_movemask_ps(_mm256_and_ps(inside, inside_depth));
                                alignas(32) float a_r[8];
                                alignas(32) float a_g[8];
                                alignas(32) float a_b[8];
                                _mm256_store_ps(a_r, a_vr);
                                _mm256_store_ps(a_g, a_vg);
                                _mm256_store_ps(a_b, a_vb);
                                //extract the depth data from the SIMD register
                                alignas(32) float depth_8f[8];
                                _mm256_store_ps(depth_8f, depth_v);
                                unsigned char r, g, b;
                                float depth_end;
                                //可以考虑顺序使用八个r,八个g,八个b

                                while (realmask) {
                                    int i = _tzcnt_u32(realmask);
                                    realmask &= (realmask - 1);
                                    r = a_r[i] * 255;
                                    g = a_g[i] * 255;
                                    b = a_b[i] * 255;
                                    depth_end = depth_8f[i];

                                    renderer.canvas.draw(x + i, y, r, g, b);
                                    renderer.zbuffer(x + i, y) = depth_end;


                                }
                                /*
                                if (realmask == 0xFF) {
                                    alignas(32) float a_r[8];
                                    alignas(32) float a_g[8];
                                    alignas(32) float a_b[8];
                                   
                                    a_vr = _mm256_mul_ps(twofivefive, a_vr);
									a_vg = _mm256_mul_ps(twofivefive, a_vg);
									a_vb = _mm256_mul_ps(twofivefive, a_vb);


                                    _mm256_store_ps(a_r, a_vr);
                                    _mm256_store_ps(a_g, a_vg);
                                    _mm256_store_ps(a_b, a_vb);

                                    
									renderer.canvas.draw(x, y, a_r, a_g, a_b);
                                    alignas(32) float depth_8f[8];
                                    _mm256_store_ps(depth_8f, depth_v);
                                    for(int i=0;i<8;i++)
                                    {
                                       
                                        renderer.zbuffer(x + i, y) = depth_8f[i];
									}
                                }
                                else {

                                    alignas(32) float a_r[8];
                                    alignas(32) float a_g[8];
                                    alignas(32) float a_b[8];
                                    _mm256_store_ps(a_r, a_vr);
                                    _mm256_store_ps(a_g, a_vg);
                                    _mm256_store_ps(a_b, a_vb);
                                    //extract the depth data from the SIMD register
                                    alignas(32) float depth_8f[8];
                                    _mm256_store_ps(depth_8f, depth_v);
                                    unsigned char r, g, b;
                                    float depth_end;
                                    //可以考虑顺序使用八个r,八个g,八个b
                                  
                                    while (realmask) {
                                        int i = _tzcnt_u32(realmask);
                                        realmask &= (realmask - 1);
                                        r = a_r[i] * 255;
                                        g = a_g[i] * 255;
                                        b = a_b[i] * 255;
                                        depth_end = depth_8f[i];

                                        renderer.canvas.draw(x + i, y, r, g, b);
                                        renderer.zbuffer(x + i, y) = depth_end;


                                    }
                                }*/
                              





                            }
                        
                    
                
					w0_vec = _mm256_add_ps(w0_vec, w0_stepx_v);
					w1_vec = _mm256_add_ps(w1_vec, w1_stepx_v);
					w2_vec = _mm256_add_ps(w2_vec, w2_stepx_v);
                   
                    //++dcdx
                    col_v_r = _mm256_add_ps(col_v_r, dcdx_block_vr);
                    col_v_g = _mm256_add_ps(col_v_g, dcdx_block_vg);
                    col_v_b = _mm256_add_ps(col_v_b, dcdx_block_vb);
                    //++dzdx
                    depth_v = _mm256_add_ps(depth_v, dzdx_block_v);
                    //++dndx
                    normal_v_x = _mm256_add_ps(normal_v_x, dndx_block_vx);
                    normal_v_y = _mm256_add_ps(normal_v_y, dndx_block_vy);
                    normal_v_z = _mm256_add_ps(normal_v_z, dndx_block_vz);
                   // z += dzdx * 8;
                   // col = col + dcdx * 8;
                   // nor = nor + dndx * 8;
                    
                    if (maxX - x < 8) {
                        __m256 m0 = _mm256_cmp_ps(w0_vec, zero, _CMP_GE_OQ);
                        __m256 m1 = _mm256_cmp_ps(w1_vec, zero, _CMP_GE_OQ);
                        __m256 m2 = _mm256_cmp_ps(w2_vec, zero, _CMP_GE_OQ);
                        __m256 inside = _mm256_and_ps(_mm256_and_ps(m0, m1), m2);
                        int mask = _mm256_movemask_ps(inside);
                        if (mask == 0) {
                            break;
                        }
                       
                        __m256 zbuffer_v = _mm256_loadu_ps(&renderer.zbuffer(x, y));
                        __m256 m0_depth = _mm256_cmp_ps(depth_v, zeor_dot_oneoneone, _CMP_GE_OQ);
                        __m256 m1_depth = _mm256_cmp_ps(zbuffer_v, depth_v, _CMP_GE_OQ);
                        __m256 inside_depth = _mm256_and_ps(m0_depth, m1_depth);
                        int mask_depth = _mm256_movemask_ps(inside_depth);

                        if (mask_depth == 0) {
                            break;
                        }

                        int realmask = _mm256_movemask_ps(_mm256_and_ps(inside, inside_depth));

                        
                           

                            col_v_r = _mm256_min_ps(col_v_r, one);
                          
                            __m256 normal_length = _mm256_sqrt_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(normal_v_x, normal_v_x), _mm256_mul_ps(normal_v_y, normal_v_y)), _mm256_mul_ps(normal_v_z, normal_v_z)));
                            normal_v_x = _mm256_div_ps(normal_v_x, normal_length);
                            normal_v_y = _mm256_div_ps(normal_v_y, normal_length);
                            normal_v_z = _mm256_div_ps(normal_v_z, normal_length);
                            
                            __m256 dot_v = _mm256_max_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(normal_v_x, light_omega_vx), _mm256_mul_ps(normal_v_y, light_omega_vy)), _mm256_mul_ps(normal_v_z, light_omega_vz)), zero);

                          
                            __m256 a_vr = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(col_v_r, light_kd), _mm256_mul_ps(light_l_vr, dot_v)), _mm256_mul_ps(light_ambient_vr, light_ka));
                            __m256 a_vg = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(col_v_g, light_kd), _mm256_mul_ps(light_l_vg, dot_v)), _mm256_mul_ps(light_ambient_vg, light_ka));
                            __m256 a_vb = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(col_v_b, light_kd), _mm256_mul_ps(light_l_vb, dot_v)), _mm256_mul_ps(light_ambient_vb, light_ka));


                            alignas(32) float a_r[8];
                            alignas(32) float a_g[8];
                            alignas(32) float a_b[8];
                            _mm256_store_ps(a_r, a_vr);
                            _mm256_store_ps(a_g, a_vg);
                            _mm256_store_ps(a_b, a_vb);
                          
                            alignas(32) float depth_8f[8];
                            _mm256_store_ps(depth_8f, depth_v);
                            unsigned char r, g, b;
                            float depth_end;
                          
                            for (int i = 0; i < maxX - x;i++) {
                               
                                r = a_r[i] * 255;
                                g = a_g[i] * 255;
                                b = a_b[i] * 255;
                                depth_end = depth_8f[i];

                                renderer.canvas.draw(x + i, y, r, g, b);
                                renderer.zbuffer(x + i, y) = depth_end;


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
    auto EdgeFunction(float xa, float ya, float xb, float yb, float x, float y) {
        return (x - xa) * (yb - ya) - (y - ya) * (xb - xa);
    }
   

    void draw_with_depth_vector(Renderer& renderer, Light& L, float ka, float kd) {//not good


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
        float w2_stepx = e2.x, w2_stepy = e2.y;
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

        alignas(32) float depth_first[8] = { z_row + dzdx * 0, z_row + dzdx * 1 , z_row + dzdx * 2 , z_row + dzdx * 3 , z_row + dzdx * 4 , z_row + dzdx * 5 , z_row + dzdx * 6 , z_row + dzdx * 7 };
        __m256 z_orw_v = _mm256_load_ps(depth_first);
        __m256 depth_blockX_v = _mm256_set1_ps(dzdx * 8);
        __m256 depth_stepY_v = _mm256_set1_ps(dzdy);
        __m256 depth_realuse_v = z_orw_v;

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

            depth_realuse_v = z_orw_v;
            colour col = c_row;
            //float z = z_row;
            vec4 nor = n_row;

            for (int x = minX; x <= maxX; x += 8) {


                __m256 m0 = _mm256_cmp_ps(w0_vec, zero, _CMP_GE_OQ);
                __m256 m1 = _mm256_cmp_ps(w1_vec, zero, _CMP_GE_OQ);
                __m256 m2 = _mm256_cmp_ps(w2_vec, zero, _CMP_GE_OQ);
                __m256 inside = _mm256_and_ps(_mm256_and_ps(m0, m1), m2);
                int mask = _mm256_movemask_ps(inside);
                if (mask == 0) {
                    w0_vec = _mm256_add_ps(w0_vec, w0_stepx_v);
                    w1_vec = _mm256_add_ps(w1_vec, w1_stepx_v);
                    w2_vec = _mm256_add_ps(w2_vec, w2_stepx_v);
                    depth_realuse_v = _mm256_add_ps(depth_realuse_v, depth_blockX_v);
                    // z += dzdx * 8;
                    col = col + dcdx * 8;
                    nor = nor + dndx * 8;
                    continue;
                }

                __m256 depth_cmp = _mm256_cmp_ps(_mm256_load_ps(&renderer.zbuffer(x, y)), depth_realuse_v, _CMP_GT_OQ);
                int depth_mask = _mm256_movemask_ps(depth_cmp);
                if (depth_mask == 0) {
                    w0_vec = _mm256_add_ps(w0_vec, w0_stepx_v);
                    w1_vec = _mm256_add_ps(w1_vec, w1_stepx_v);
                    w2_vec = _mm256_add_ps(w2_vec, w2_stepx_v);
                    depth_realuse_v = _mm256_add_ps(depth_realuse_v, depth_blockX_v);
                    // z += dzdx * 8;
                    col = col + dcdx * 8;
                    nor = nor + dndx * 8;
                    continue;
                }
                float z[8];
                _mm256_store_ps(z, depth_realuse_v);
                while (mask) {
                    int i = _tzcnt_u32(mask);
                    mask &= (mask - 1);






                    float depth = z[i];


                    // if (renderer.zbuffer(x+i, y) > depth && depth > 0.001f) {

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
                    // }

                }
                w0_vec = _mm256_add_ps(w0_vec, w0_stepx_v);
                w1_vec = _mm256_add_ps(w1_vec, w1_stepx_v);
                w2_vec = _mm256_add_ps(w2_vec, w2_stepx_v);
                depth_realuse_v = _mm256_add_ps(depth_realuse_v, depth_blockX_v);
                col = col + dcdx * 8;
                nor = nor + dndx * 8;

                if (maxX - x < 8) {
                    _mm256_store_ps(w0, w0_vec);
                    _mm256_store_ps(w1, w1_vec);
                    _mm256_store_ps(w2, w2_vec);
                    for (int i = 0; i <= maxX - x; i++) {
                        if (w0[i] >= 0 && w1[i] >= 0 && w2[i] >= 0) {
                            float depth = z[i];


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
            z_orw_v = _mm256_add_ps(z_orw_v, depth_stepY_v);
            //z_row += dzdy;
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
