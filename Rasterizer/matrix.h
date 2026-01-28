#pragma once

#include <iostream>
#include <vector>
#include "vec4.h"
#include <immintrin.h>
using namespace std;
// Matrix class for 4x4 transformation matrices
class matrix {
public:
    union {
        alignas(16) float m[4][4]; // 2D array representation of the matrix
        alignas(16) float a[16];   // 1D array representation of the matrix for linear access
    };


    // Default constructor initializes the matrix as an identity matrix
    matrix() {
        identity();
    }

    // Access matrix elements by row and column
    float& operator()(unsigned int row, unsigned int col) { return m[row][col]; }

    // Display the matrix elements in a readable format
    void display() {
        for (unsigned int i = 0; i < 4; i++) {
            for (unsigned int j = 0; j < 4; j++)
                std::cout << m[i][j] << '\t';
            std::cout << std::endl;
        }
    }

    // Multiply the matrix by a 4D vector
    // Input Variables:
    // - v: vec4 object to multiply with the matrix
    // Returns the resulting transformed vec4
    vec4 operator * (const vec4& v) const {
        vec4 result;
		__m128 value = _mm_loadu_ps(v.v);
		__m128 row0 = _mm_loadu_ps(&a[0]);
		__m128 row1 = _mm_loadu_ps(&a[4]);
		__m128 row2 = _mm_loadu_ps(&a[8]);
		__m128 row3 = _mm_loadu_ps(&a[12]);
		result.v[0] = _mm_cvtss_f32(_mm_dp_ps(row0, value, 0xF1));//dot product 
		result.v[1] = _mm_cvtss_f32(_mm_dp_ps(row1, value, 0xF1));
		result.v[2] = _mm_cvtss_f32(_mm_dp_ps(row2, value, 0xF1));
		result.v[3] = _mm_cvtss_f32(_mm_dp_ps(row3, value, 0xF1));
       
        return result;
    }

    // Multiply the matrix by another matrix
    // Input Variables:
    // - mx: Another matrix to multiply with
    // Returns the resulting matrix
    matrix operator * (const matrix& mx) const {
        matrix ret;
        __m128 rowVec1 = _mm_loadu_ps(&a[0]);
        __m128 rowVec2 = _mm_loadu_ps(&a[4]);
        __m128 rowVec3 = _mm_loadu_ps(&a[8]);
        __m128 rowVec4 = _mm_loadu_ps(&a[12]);

       
            for (int col = 0; col < 4; ++col) {
				
				__m128 colVec = _mm_set_ps(mx.a[3 * 4 + col], mx.a[2 * 4 + col], mx.a[1 * 4 + col], mx.a[0 * 4 + col]);

				ret.a[0 *4+ col] = _mm_cvtss_f32(_mm_dp_ps(rowVec1, colVec, 0xF1));
                ret.a[1 * 4 + col] = _mm_cvtss_f32(_mm_dp_ps(rowVec2, colVec, 0xF1));
				ret.a[2 * 4 + col] = _mm_cvtss_f32(_mm_dp_ps(rowVec3, colVec, 0xF1));
                ret.a[3 * 4 + col] = _mm_cvtss_f32(_mm_dp_ps(rowVec4, colVec, 0xF1));
                    
            }
        
        return ret;
    }

    // Create a perspective projection matrix
    // Input Variables:
    // - fov: Field of view in radians
    // - aspect: Aspect ratio of the viewport
    // - n: Near clipping plane
    // - f: Far clipping plane
    // Returns the perspective matrix
    static matrix makePerspective(float fov, float aspect, float n, float f) {
        matrix m;
        m.zero();
        float tanHalfFov = std::tan(fov / 2.0f);

        m.a[0] = 1.0f / (aspect * tanHalfFov);
        m.a[5] = 1.0f / tanHalfFov;
        m.a[10] = -f / (f - n);
        m.a[11] = -(f * n) / (f - n);
        m.a[14] = -1.0f;
        return m;
    }

    // Create a translation matrix
    // Input Variables:
    // - tx, ty, tz: Translation amounts along the X, Y, and Z axes
    // Returns the translation matrix
    static matrix makeTranslation(float tx, float ty, float tz) {
        matrix m;
        m.identity();
        m.a[3] = tx;
        m.a[7] = ty;
        m.a[11] = tz;
        return m;
    }

    // Create a rotation matrix around the Z-axis
    // Input Variables:
    // - aRad: Rotation angle in radians
    // Returns the rotation matrix
    static matrix makeRotateZ(float aRad) {
        matrix m;
        m.identity();
        m.a[0] = std::cos(aRad);
        m.a[1] = -std::sin(aRad);
        m.a[4] = std::sin(aRad);
        m.a[5] = std::cos(aRad);
        return m;
    }

    // Create a rotation matrix around the X-axis
    // Input Variables:
    // - aRad: Rotation angle in radians
    // Returns the rotation matrix
    static matrix makeRotateX(float aRad) {
        matrix m;
        m.identity();
        m.a[5] = std::cos(aRad);
        m.a[6] = -std::sin(aRad);
        m.a[9] = std::sin(aRad);
        m.a[10] = std::cos(aRad);
        return m;
    }

    // Create a rotation matrix around the Y-axis
    // Input Variables:
    // - aRad: Rotation angle in radians
    // Returns the rotation matrix
    static matrix makeRotateY(float aRad) {
        matrix m;
        m.identity();
        m.a[0] = std::cos(aRad);
        m.a[2] = std::sin(aRad);
        m.a[8] = -std::sin(aRad);
        m.a[10] = std::cos(aRad);
        return m;
    }

    // Create a composite rotation matrix from X, Y, and Z rotations
    // Input Variables:
    // - x, y, z: Rotation angles in radians around each axis
    // Returns the composite rotation matrix
    static matrix makeRotateXYZ(float x, float y, float z) {
        float sx = sinf(x), cx = cosf(x);
    float sy = sinf(y), cy = cosf(y);
    float sz = sinf(z), cz = cosf(z);

    matrix m;

    // R = Rz * Ry * Rx （常用约定）
    m.a[0]  =  cy * cz;
    m.a[1]  =  cy * sz;
    m.a[2]  = -sy;
    m.a[3]  =  0.0f;

    m.a[4]  =  sx * sy * cz - cx * sz;
    m.a[5]  =  sx * sy * sz + cx * cz;
    m.a[6]  =  sx * cy;
    m.a[7]  =  0.0f;

    m.a[8]  =  cx * sy * cz + sx * sz;
    m.a[9]  =  cx * sy * sz - sx * cz;
    m.a[10] =  cx * cy;
    m.a[11] =  0.0f;

    m.a[12] =  0.0f;
    m.a[13] =  0.0f;
    m.a[14] =  0.0f;
    m.a[15] =  1.0f;

    return m;
    }

    // Create a scaling matrix
    // Input Variables:
    // - s: Scaling factor
    // Returns the scaling matrix
    static matrix makeScale(float s) {
        matrix m;
        s = std::max(s, 0.01f); // Ensure scaling factor is not too small
        m.identity();
        m.a[0] = s;
        m.a[5] = s;
        m.a[10] = s;
        return m;
    }

    // Create an identity matrix
    // Returns an identity matrix
    static matrix makeIdentity() {
        matrix m;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                m.m[i][j] = (i == j) ? 1.0f : 0.0f;
            }
        }
        return m;
    }

private:
    // Set all elements of the matrix to 0
    void zero() {
        for (unsigned int i = 0; i < 16; i++)
            a[i] = 0.f;
    }

    // Set the matrix as an identity matrix
    void identity() {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                m[i][j] = (i == j) ? 1.0f : 0.0f;
            }
        }
    }
};


