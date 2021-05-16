#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Task 7: Implement this

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  // fill all 0 sub levels
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {
    fill_mip_level(tex, i);
  }
}


int prev_mip_loc(int mip_loc, int boundary) {
  return 2 * mip_loc;
}


int prev_mip_loc_boundary_check(int mip_loc, int boundary) {
  return min(2 * mip_loc, boundary - 1);
}

void Sampler2DImp::fill_mip_level(Texture& tex, int level) {
  if (level < 1) {
    return;
  }

  auto prev_fct = prev_mip_loc;

  if (level == 1) {
    prev_fct = prev_mip_loc_boundary_check;
  }

  MipLevel& mip = tex.mipmap[level];
  MipLevel& prev_mip = tex.mipmap[level - 1];

  int t, t_prev, i_prev, j_prev;
  int r = 0, g = 0, b = 0, a = 0;

  for (size_t i = 0; i < mip.width; ++i) {
    for (size_t j = 0; j < mip.height; ++j) {

      t = 4 * (i + j * mip.width);

      i_prev = prev_fct(i, prev_mip.width);
      j_prev = prev_fct(j, prev_mip.height);

      t_prev = 4 * (i_prev + j_prev * prev_mip.width);

      r = 0, g = 0, b = 0, a = 0;

      r += prev_mip.texels[t_prev    ];
      g += prev_mip.texels[t_prev + 1];
      b += prev_mip.texels[t_prev + 2];
      a += prev_mip.texels[t_prev + 3];

      t_prev += 4; // i++

      r += prev_mip.texels[t_prev    ];
      g += prev_mip.texels[t_prev + 1];
      b += prev_mip.texels[t_prev + 2];
      a += prev_mip.texels[t_prev + 3];

      t_prev += 4 * prev_mip.width; // j++

      r += prev_mip.texels[t_prev    ];
      g += prev_mip.texels[t_prev + 1];
      b += prev_mip.texels[t_prev + 2];
      a += prev_mip.texels[t_prev + 3];

      t_prev -= 4; // i--

      r += prev_mip.texels[t_prev    ];
      g += prev_mip.texels[t_prev + 1];
      b += prev_mip.texels[t_prev + 2];
      a += prev_mip.texels[t_prev + 3];

      mip.texels[t    ] = (unsigned char)round(r / 4.0f);
      mip.texels[t + 1] = (unsigned char)round(g / 4.0f);
      mip.texels[t + 2] = (unsigned char)round(b / 4.0f);
      mip.texels[t + 3] = (unsigned char)round(a / 4.0f);
    }
  }
}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  // Task 6: Implement nearest neighbour interpolation

  // mip level
  if (level >= tex.mipmap.size()) {
    // return magenta for invalid level
    return Color(1, 0, 1, 1);
  }

  auto& mip = tex.mipmap[level];

  // map u,v [0,1] to nearest image pixel
  int tx = (int)floor(u * mip.width);
  int ty = (int)floor(v * mip.height);

  int t = 4 * (tx + ty * mip.width);

  return Color(
    mip.texels[t    ] / 255.0f,
    mip.texels[t + 1] / 255.0f,
    mip.texels[t + 2] / 255.0f,
    mip.texels[t + 3] / 255.0f
  );
}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task 6: Implement bilinear filtering

  // mip level
  if (level >= tex.mipmap.size()) {
    // return magenta for invalid level
    return Color(1, 0, 1, 1);
  }

  auto& mip = tex.mipmap[level];

  // map to image
  u *= mip.width;
  v *= mip.height;

  // get closest image loc
  int i = (int)floor(u - 0.5f);
  int j = (int)floor(v - 0.5f);

  if (i < 0) {
    i = 0;
  }

  if (j < 0) {
    j = 0;
  }

  int at = 4 * (i + j * mip.width);

  auto c00 = Color(
    mip.texels[at    ] / 255.0f,
    mip.texels[at + 1] / 255.0f,
    mip.texels[at + 2] / 255.0f,
    mip.texels[at + 3] / 255.0f
  );

  if (i < mip.width - 1) {
    at += 4; // tx++
  }

  auto c10 = Color(
    mip.texels[at    ] / 255.0f,
    mip.texels[at + 1] / 255.0f,
    mip.texels[at + 2] / 255.0f,
    mip.texels[at + 3] / 255.0f
  );

  if (j < mip.height - 1) {
    at += 4 * mip.width; // j++
  }
  
  auto c11 = Color(
    mip.texels[at    ] / 255.0f,
    mip.texels[at + 1] / 255.0f,
    mip.texels[at + 2] / 255.0f,
    mip.texels[at + 3] / 255.0f
  );

  if (i < mip.width - 1) {
    at -= 4; // i--
  }

  auto c01 = Color(
    mip.texels[at    ] / 255.0f,
    mip.texels[at + 1] / 255.0f,
    mip.texels[at + 2] / 255.0f,
    mip.texels[at + 3] / 255.0f
  );
  
  float s = u - (i + 0.5);
  float t = v - (j + 0.5);

  return (1.0 - t) * ((1 - s) * c00 + s * c10) + t * ((1.0 - s) * c01 + s * c11);
}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Task 7: Implement trilinear filtering

  float d = log2f(max(u_scale, v_scale));

  // find 2 levels to interpolate
  int level = (int)floor(d);

  if (level < 0) {
    return sample_bilinear(tex, u, v, 0);
  }
  else if (level >= tex.mipmap.size() - 1) {
    return sample_bilinear(tex, u, v, tex.mipmap.size() - 1);
  }

  auto c0 = sample_bilinear(tex, u, v, level);
  auto c1 = sample_bilinear(tex, u, v, level + 1);

  float w = d - level;

  return (1.0 - w) * c0 + w * c1;
}

} // namespace CMU462
