#pragma once
#include <thread>
#include <vector>
#include <functional>
#include <memory>
#include "colour.h"
#include "renderer.h"
#include "light.h"
#include <mutex>
#include <atomic>
#include <queue>
#include <chrono>
struct alignas(64) TileFlag {
	std::atomic<bool> free;
};
//struct for data in every single tile
struct alignas(32) TileWork {
    Renderer* renderer;
    Light light;

    int minX, maxX, minY, maxY, ydifferent;

    float row_w0, row_w1, row_w2;
    __m256 w0_step_vx_256, w1_step_vx_256, w2_step_vx_256;

    colour c_row;
    __m256 dcdx_v_r, dcdx_v_g, dcdx_v_b;

    float z_row;
    __m256 dzdx_v_r;

    vec4 n_row;
    __m256 dndx_v_x, dndx_v_y, dndx_v_z;

    __m256 w0_stepx_v, w1_stepx_v, w2_stepx_v;

    __m256 dcdx_block_vr, dcdx_block_vg, dcdx_block_vb;
    __m256 dzdx_block_v;
    __m256 dndx_block_vx, dndx_block_vy, dndx_block_vz;

    float w0_stepy, w1_stepy, w2_stepy;
    float dzdy;
    colour dcdy;
    vec4 dndy;
    float ka, kd;
};


class SPSCQueue {
public:
	
	//task queue for one tile
    std::queue<TileWork> taskQueue;
 
    bool try_push(const TileWork& v) {
       
        taskQueue.push(v);
        return true;
    }

   

    bool try_pop(TileWork& out) {
     ;
        if (taskQueue.empty()) {
            return false; // Queue is empty,work end 
		}
        out = taskQueue.front();
        taskQueue.pop();
		return true;
    }
	
	
};


class MultilThreadControl
{
	public:
        std::vector<double> tile_draw_number;
        std::atomic<int> stop_flag = 0;
		std::mutex tex;
		std::atomic<bool> one_done = false;
		std::atomic<int>active_workers=0;
        std::map<int,int>massion_owner;
        
       
		SPSCQueue* tiles = new SPSCQueue[32];
        //vector<SPSCQueue> tiles;
        //SPSCQueue tiles[32];
		int tile_count = 8;
		//std::vector<TileFlag>is_freetile;
		std::atomic<bool> produce_done = false;
	static const int MAX_THREADS = 10;
	 int numThreads=0;
	std::vector<std::jthread> scv;
	MultilThreadControl() {
		
	}
	void start(int n=10) {
        numThreads=n;
		scv.reserve(n);
        tile_draw_number.reserve(n);
		for (int i =1; i < n; i++)
		{
			massion_owner[i] = i;
            tile_draw_number.emplace_back(0.f); // Records each worker thread's per-frame execution time
            scv.emplace_back(&MultilThreadControl::worker, this, i);
			
		}
        
	}
  
	void setTileCount(int n) {
        // Sets the number of screen tiles
        tile_count = n;
	}
	
    inline void tile_draw(Renderer& renderer, int minX, int maxX, int minY, int maxY, int ydifferent, float row_w0, float row_w1, float row_w2, __m256 w0_step_vx_256
        , __m256 w1_step_vx_256, __m256 w2_step_vx_256, colour c_row, __m256 dcdx_v_r, __m256 dcdx_v_g, __m256 dcdx_v_b, float z_row, __m256 dzdx_v_r,
        vec4 n_row, __m256 dndx_v_x, __m256 dndx_v_y, __m256 dndx_v_z, __m256 w0_stepx_v, __m256 w1_stepx_v, __m256 w2_stepx_v,
        __m256 dcdx_block_vr, __m256 dcdx_block_vg, __m256 dcdx_block_vb, __m256 dzdx_block_v, __m256 dndx_block_vx, __m256 dndx_block_vy, __m256 dndx_block_vz,
        float w0_stepy, float w1_stepy, float w2_stepy, float dzdy, colour dcdy, vec4 dndy, float ka, float kd, Light L)
    {
        __m256 one = _mm256_set1_ps(1.0f);
        __m256 zeor_dot_oneoneone = _mm256_set1_ps(0.001f);
        __m256 twofivefive = _mm256_set1_ps(255.0f);
        __m256 zero = _mm256_setzero_ps();

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


        row_w0 += w0_stepy * ydifferent;
        row_w1 += w1_stepy * ydifferent;
        row_w2 += w2_stepy * ydifferent;
        z_row += dzdy * ydifferent;
        c_row = c_row + dcdy * ydifferent;
        n_row = n_row + dndy * ydifferent;
        // To be implemented: tile-based rasterization

        for (int y = minY; y <= maxY; ++y) {
            if (maxY > 1000) {
                int xsc = 1;
            }
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

                __m256 zbuffer_v = _mm256_loadu_ps(&Renderer::instance().zbuffer(x, y));
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


                        Renderer::instance().canvas.draw(x, y, a_r, a_g, a_b);

                        _mm256_store_ps(&Renderer::instance().zbuffer(x, y), depth_v);

                    }
                    else {

                        alignas(32) float a_r[8];
                        alignas(32) float a_g[8];
                        alignas(32) float a_b[8];
                        a_vr = _mm256_mul_ps(twofivefive, a_vr);
                        a_vg = _mm256_mul_ps(twofivefive, a_vg);
                        a_vb = _mm256_mul_ps(twofivefive, a_vb);
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
                            r = a_r[i];
                            g = a_g[i];
                            b = a_b[i];
                            depth_end = depth_8f[i];

                            Renderer::instance().canvas.draw(x + i, y, r, g, b);
                            Renderer::instance().zbuffer(x + i, y) = depth_end;


                        }
                    }






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

                    __m256 zbuffer_v = _mm256_loadu_ps(&Renderer::instance().zbuffer(x, y));
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

                    for (int i = 0; i < maxX - x; i++) {

                        r = a_r[i] * 255;
                        g = a_g[i] * 255;
                        b = a_b[i] * 255;
                        depth_end = depth_8f[i];

                        Renderer::instance().canvas.draw(x + i, y, r, g, b);
                        Renderer::instance().zbuffer(x + i, y) = depth_end;

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
	

private:
	void worker(int tid) {
		TileWork mission;
		
        int my_epoch = 0;
		while (true)
		{
            //int my_epoch = stop_flag.load(std::memory_order_acquire);
            // Wait for the main thread to signal the start of a new frame 
            stop_flag.wait(my_epoch);
            my_epoch = stop_flag.load(std::memory_order_acquire);
           
            auto star2 = std::chrono::high_resolution_clock::now();//count time
            // The first dequeued item is used to obtain the tile's minY/maxY range for clearing buffers
            if (tiles[tid].try_pop(mission)) {
                Renderer::instance().zbuffer.tile_clear(mission.minY, mission.maxY);
                Renderer::instance().canvas.tile_clear(mission.minY, mission.maxY);
            }
          

            
                // Consume and rasterize all tasks for the tile owned by this thread
                while (tiles[tid].try_pop(mission)) {
                    // Renderer has been refactored into a global singleton.
                    // The pointer will still be passed, but in fact it is no longer being used. I just too lazy to change it.
                    tile_draw(*mission.renderer, mission.minX, mission.maxX, mission.minY, mission.maxY, mission.ydifferent,
                        mission.row_w0, mission.row_w1, mission.row_w2, mission.w0_step_vx_256,
                        mission.w1_step_vx_256, mission.w2_step_vx_256, mission.c_row, mission.dcdx_v_r, mission.dcdx_v_g, mission.dcdx_v_b,
                        mission.z_row, mission.dzdx_v_r,
                        mission.n_row, mission.dndx_v_x, mission.dndx_v_y, mission.dndx_v_z, mission.w0_stepx_v, mission.w1_stepx_v, mission.w2_stepx_v,
                        mission.dcdx_block_vr, mission.dcdx_block_vg, mission.dcdx_block_vb, mission.dzdx_block_v,
                        mission.dndx_block_vx, mission.dndx_block_vy, mission.dndx_block_vz,
                        mission.w0_stepy, mission.w1_stepy, mission.w2_stepy, mission.dzdy, mission.dcdy, mission.dndy,
                        mission.ka, mission.kd, mission.light);
                    
				}
                auto end2 = std::chrono::high_resolution_clock::now();// End timing
                tile_draw_number[tid] = std::chrono::duration<double, std::milli>(end2 - star2).count();
                {
					
                    
                  
                    // active_workers is an atomic counter:
                    // when active_workers reaches 0, it means all worker threads have finished this frame
                    active_workers--;
					
                  

                   
                }
				
            
                            
			
		}
	}



   




	
};

