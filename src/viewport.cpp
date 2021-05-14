#include "viewport.h"

#include "CMU462.h"

namespace CMU462 {

void ViewportImp::set_viewbox( float centerX, float centerY, float vspan ) {

  // Task 5 (part 2): 
  // Set svg coordinate to normalized device coordinate transformation. Your input
  // arguments are defined as normalized SVG canvas coordinates.
  this->centerX = centerX;
  this->centerY = centerY;
  this->vspan = vspan; 

  // m takes (centerX, centerY) to (0.5, 0.5)
  // m takes (centerX, centerY + vspan) to (0.5, 1.0)
  // m preserves aspect ratio

  /*
    s  0  tx
    0  s  ty
    0  0  1
  */
  Matrix3x3 m = Matrix3x3::identity();

  float s = 0.5f / vspan;
  float tx = 0.5f - s * centerX;
  float ty = 0.5f - s * centerY;

  m(0, 0) = s;
  m(1, 1) = s;
  
  m(0, 2) = tx;
  m(1, 2) = ty;

  set_svg_2_norm(m);
}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->centerX -= dx;
  this->centerY -= dy;
  this->vspan *= scale;
  set_viewbox( centerX, centerY, vspan );
}

} // namespace CMU462
