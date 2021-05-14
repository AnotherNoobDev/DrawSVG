#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <utility>
#include <cassert>

#include "triangulation.h"

const float VERTICAL_LINE_EPSILON = 0.001;
const float FLOAT_POINT_EPSILON = 0.001;

using namespace std;

namespace CMU462 {


// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg( SVG& svg ) {

  // set top level transformation
  transformation = svg_2_screen;

  // draw all elements
  for ( size_t i = 0; i < svg.elements.size(); ++i ) {
    draw_element(svg.elements[i]);
  }

  // draw canvas outline
  Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y--;
  Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y++;

  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  resolve();

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;

}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;

}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

  // Task 5 (part 1):
  // Modify this to implement the transformation stack

  switch(element->type) {
    case POINT:
      draw_point(static_cast<Point&>(*element));
      break;
    case LINE:
      draw_line(static_cast<Line&>(*element));
      break;
    case POLYLINE:
      draw_polyline(static_cast<Polyline&>(*element));
      break;
    case RECT:
      draw_rect(static_cast<Rect&>(*element));
      break;
    case POLYGON:
      draw_polygon(static_cast<Polygon&>(*element));
      break;
    case ELLIPSE:
      draw_ellipse(static_cast<Ellipse&>(*element));
      break;
    case IMAGE:
      draw_image(static_cast<Image&>(*element));
      break;
    case GROUP:
      draw_group(static_cast<Group&>(*element));
      break;
    default:
      break;
  }

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Extra credit 

}

void SoftwareRendererImp::draw_image( Image& image ) {

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

int SoftwareRendererImp::closest_cell(float x) {
  /*
  int res = (int)round(x);

  // for 0.5, round down
  if (abs(wu_fpart(x) - 0.5) < FLOAT_POINT_EPSILON) {
    res = (int)floor(x);
  }

  return res;
  */

  /*
  if (wu_fpart(x) > 0.5) {
    return (int)ceil(x);
  }
  else {
    return floor(x);
  }
  */

  return floor(x);
}

SoftwareRendererImp::SamplingRange SoftwareRendererImp::get_sampling_range(float x0, float x1, float upper_bound) {
  assert(x0 < x1 + FLOAT_POINT_EPSILON);

  float sample_dist = 1.0;
  float start, end;

  // assuming 1 sample per pixel for now
  float x0_fpart = wu_fpart(x0);
  float x0_ipart = wu_ipart(x0);

  start = x0_ipart + 0.5;

  // check if x0 is to the right of the sample
  if (x0_fpart > 0.5 + FLOAT_POINT_EPSILON) {
    start += sample_dist; //next sample
  }

  if (start < 0.5) {
    start = 0.5;
  }

  float x1_fpart = wu_fpart(x1);
  float x1_ipart = wu_ipart(x1);

  end = x1_ipart + 0.5;
  
  // check if x1 is to the left of the sample
  if (x1_fpart < 0.5 - FLOAT_POINT_EPSILON) {
    end -= sample_dist; //prev sample
  }

  if (end > upper_bound - 0.5) {
    end = upper_bound - 0.5;
  }

  end += sample_dist;

  return SamplingRange(start, end, sample_dist);
}

void SoftwareRendererImp::fill_point(int sx, int sy, Color color) {
  // fill sample - NOT doing alpha blending!
  render_target[4 * (sx + sy * target_w)] = (uint8_t)(color.r * 255);
  render_target[4 * (sx + sy * target_w) + 1] = (uint8_t)(color.g * 255);
  render_target[4 * (sx + sy * target_w) + 2] = (uint8_t)(color.b * 255);
  render_target[4 * (sx + sy * target_w) + 3] = (uint8_t)(color.a * 255);
}

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int) floor(x);
  int sy = (int) floor(y);

  // check bounds
  if ( sx < 0 || sx >= target_w ) return;
  if ( sy < 0 || sy >= target_h ) return;

  fill_point(sx, sy, color);
}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {

  // Task 2: 
  // Implement line rasterization

  // TODO(optional) line width (strokeWidth)
  // idea: 
  // rasterize aliased line of width - 2; 
  // use wu for outer edges

  // TODO don't think it matches reference that well, check at the end
  // check images 2,3 (basic)
  
  bool steep = abs(y1 - y0) > abs(x1 - x0);
  
  if (steep) {
    swap(x0, y0);
    swap(x1, y1);
  }
      
  if (x0 > x1) {
    swap(x0, x1);
    swap(y0, y1);
  }

  float dx = x1 - x0;
  float dy = y1 - y0;
  float gradient;

  if (dx < VERTICAL_LINE_EPSILON) {
    gradient = 1.0;
  }
  else {
    gradient = dy / dx;
  }
  
  // handle first endpoint
  float xend = wu_round(x0);
  float yend = y0 + gradient * (xend - x0);
  float xgap = wu_rfpart(x0 + 0.5);
  float xpxl1 = xend; // this will be used in the main loop
  float ypxl1 = wu_ipart(yend);
          
  if (steep) {
    rasterize_point(ypxl1, xpxl1, (wu_rfpart(yend) * xgap) * color);
    rasterize_point(ypxl1 + 1, xpxl1, (wu_fpart(yend) * xgap) * color);
  }
  else {
    rasterize_point(xpxl1, ypxl1, (wu_rfpart(yend) * xgap) * color);
    rasterize_point(xpxl1, ypxl1 + 1, (wu_fpart(yend) * xgap) * color);
  }

  float intery = yend + gradient; // first y-intersection for the main loop
                                  
  // handle second endpoint
  xend = wu_round(x1);
  yend = y1 + gradient * (xend - x1);
  xgap = wu_fpart(x1 + 0.5);
  float xpxl2 = xend; //this will be used in the main loop
  float ypxl2 = wu_ipart(yend);

  if (steep) {
    rasterize_point(ypxl2, xpxl2, (wu_rfpart(yend) * xgap) * color);
    rasterize_point(ypxl2 + 1, xpxl2, (wu_fpart(yend) * xgap) * color);
  }
  else {
    rasterize_point(xpxl2, ypxl2, (wu_rfpart(yend) * xgap) * color);
    rasterize_point(xpxl2, ypxl2 + 1, (wu_fpart(yend) * xgap) * color);
  }
  
  // main loop
  if (steep) {
    for (float x = xpxl1 + 1; x <= xpxl2 - 1; ++x) {
      rasterize_point(wu_ipart(intery), x, wu_rfpart(intery) * color);
      rasterize_point(wu_ipart(intery) + 1, x, wu_fpart(intery) * color);
      intery = intery + gradient;
    }
  }
  else {
    for (float x = xpxl1 + 1; x <= xpxl2 - 1; ++x) {
      rasterize_point(x, wu_ipart(intery), wu_rfpart(intery) * color);
      rasterize_point(x, wu_ipart(intery) + 1, wu_fpart(intery) * color);
      intery = intery + gradient;
    }
  }
}


void SoftwareRendererImp::rasterize_triangle(float x0, float y0,
                                             float x1, float y1,
                                             float x2, float y2,
                                             Color color) {
  // Task 3: 
  // Implement triangle rasterization

  //printf("draw triangle: %f %f %f\n", x0, x1, x2);

  // sort in increasing x order

  if (x0 > x1) {
    swap(x0, x1);
    swap(y0, y1);
  }

  if (x1 > x2) {
    swap(x1, x2);
    swap(y1, y2);
  }

  if (x0 > x1) {
    swap(x0, x1);
    swap(y0, y1);
  }

  assert(x0 <= x1);
  assert(x1 <= x2);
  assert(x0 < x2);

  // slopes
  float m01 = 1.0;
  bool infSlope01 = true;
  if ((x1 - x0) > VERTICAL_LINE_EPSILON) {
    m01 = (y1 - y0) / (x1 - x0);
    infSlope01 = false;
  }

  float m12 = 1.0;
  bool infSlope02 = true;
  if ((x2 - x1) > VERTICAL_LINE_EPSILON) {
    m12 = (y2 - y1) / (x2 - x1);
    infSlope02 = false;
  }

  float m02 = (y2 - y0) / (x2 - x0);

  auto xRange = get_sampling_range(x0, x2, target_w);

  float x;
  float yBottom, yTop;


  // [x0, x1)
  for (x = xRange.start; x < x1; x += xRange.step) {
    assert(!infSlope01);
    yBottom = y0 + m02 * (x - x0);
    yTop = y0 + m01 * (x - x0);

    if (yTop < yBottom) {
      swap(yTop, yBottom);
    }

    auto yRange = get_sampling_range(yBottom, yTop, target_h);

    for (float y = yRange.start; y < yRange.stop; y += yRange.step) {
      //fill_point(closest_cell(x), closest_cell(y), Color(0.0, y / yRange.stop, 0.0));
      fill_point(closest_cell(x), closest_cell(y), color);
      //fill_point(closest_cell(x), closest_cell(y), Color(1.0,0.0,0.0));
    }
  }

  // x1
  if (abs(x1 - x) < FLOAT_POINT_EPSILON) {
    yBottom = y0 + m02 * (x1 - x0);
    yTop = y1;

    if (yTop < yBottom) {
      swap(yTop, yBottom);
    }

    auto yRange = get_sampling_range(yBottom, yTop, target_h);

    for (float y = yRange.start; y < yRange.stop; y += yRange.step) {
      fill_point(closest_cell(x), closest_cell(y), color);
      //fill_point(closest_cell(x), closest_cell(y), Color(0.0, 1.0, 0.0));
    }

    x += xRange.step;
  }

  // (x1, x2]
  for (; x < xRange.stop; x += xRange.step) {
    assert(!infSlope02);
    yBottom = y0 + m02 * (x - x0);
    yTop = y1 + m12 * (x - x1);

    if (yTop < yBottom) {
      swap(yTop, yBottom);
    }

    auto yRange = get_sampling_range(yBottom, yTop, target_h);

    for (float y = yRange.start; y < yRange.stop; y += yRange.step) {
      //fill_point(closest_cell(x), closest_cell(y), Color(y / yRange.stop, 0.0, 0.0));
      fill_point(closest_cell(x), closest_cell(y), color);
      //fill_point(closest_cell(x), closest_cell(y), Color(0.0, 0.0, 1.0));
    }
  }
}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 6: 
  // Implement image rasterization

}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {

  // Task 4: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 4".
  return;

}


} // namespace CMU462
