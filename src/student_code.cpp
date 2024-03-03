#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    // TODO Part 1.
    vector<Vector2D> next_curve;
    int num_points = points.size();
    for(int i = 0; i < num_points - 1; i++){
      Vector2D next_point = points[i] * (1 - t) + points[i+1] * t;
      next_curve.push_back(next_point);
    }
    return next_curve;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
    vector<Vector3D> next_patch;
    int num_points = points.size();
    for(int i = 0; i < num_points - 1; i++){
      Vector3D next_point = points[i] * (1 - t) + points[i+1] * t;
      next_patch.push_back(next_point);
    }
    return next_patch;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
    vector<Vector3D> next_patch;
    next_patch = evaluateStep(points, t);
    while(next_patch.size() != 1){
      next_patch = evaluateStep(next_patch, t);
    }
    return next_patch[0];
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    // TODO Part 2.
    vector<Vector3D> moving_curve;
    int num_col = controlPoints.size();
    Vector3D moving_point;
    for(int i = 0; i < num_col; i++){
      moving_point = evaluate1D(controlPoints[i], u);
      moving_curve.push_back(moving_point);
    }
    return evaluate1D(moving_curve, v);
  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.
    Vector3D nor;
    HalfedgeCIter h = this->halfedge();      // get the outgoing half-edge of the vertex
    do {
      HalfedgeCIter h_twin = h->twin(); // get the opposite half-edge
      VertexCIter v = h_twin->vertex(); // vertex is the 'source' of the half-edge, so
                                        // h->vertex() is v, whereas h_twin->vertex()
                                        // is the neighboring vertex
      h = h_twin->next();               // move to the next outgoing half-edge of the vertex
      if(!h->isBoundary())
        nor += h->face()->normal();
    } while(h != this->halfedge()); // keep going until we are back where we were
    
    return nor.unit();
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.
    auto h0 = e0->halfedge();
    auto h1 = h0->next();
    auto h2 = h1->next();
    auto h3 = h0->twin();
    auto h4 = h3->next();
    auto h5 = h4->next();

    auto v0 = h4->vertex();
    auto v1 = h1->vertex();
    auto v2 = h2->vertex();
    auto v3 = h5->vertex();

    auto e1 = h1->edge();
    auto e2 = h2->edge();
    auto e3 = h4->edge();
    auto e4 = h5->edge();

    auto f0 = h0->face();
    auto f1 = h3->face();

    if(h0->isBoundary() || h1->isBoundary() || h2->isBoundary()||
       h3->isBoundary() || h4->isBoundary() || h5->isBoundary())return e0;
    //next, twin, vertex, edge, face
    h0->setNeighbors(h5, h3, v2, e0, f0);
    h1->setNeighbors(h0, h1->twin(), v1, e1, f0);
    h2->setNeighbors(h4, h2->twin(), v2, e2, f1);
    h3->setNeighbors(h2, h0, v3, e0, f1);
    h4->setNeighbors(h3, h4->twin(), v0, e3, f1);
    h5->setNeighbors(h1, h5->twin(), v3, e4, f0);

    v0->halfedge() = h4;
    v1->halfedge() = h1;
    v2->halfedge() = h2;
    v3->halfedge() = h3;

    e0->halfedge() = h0;

    f0->halfedge() = h0;
    f1->halfedge() = h3;
    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.
    if(e0->isBoundary())
      return e0->halfedge()->vertex();
    
    auto h0 = e0->halfedge();
    auto h1 = h0->next();
    auto h2 = h1->next();
    auto h3 = h0->twin();
    auto h4 = h3->next();
    auto h5 = h4->next();

    auto v0 = h4->vertex();
    auto v1 = h1->vertex();
    auto v2 = h2->vertex();
    auto v3 = h5->vertex();

    auto e1 = h1->edge();
    auto e2 = h2->edge();
    auto e3 = h4->edge();
    auto e4 = h5->edge();

    auto f0 = h0->face();
    auto f1 = h3->face();

    auto h6 = newHalfedge(), h7 = newHalfedge(), h8 = newHalfedge(),
         h9 = newHalfedge(), h10 = newHalfedge(), h11 = newHalfedge();

    auto m = newVertex();

    auto e6 = newEdge(), e7 = newEdge(), e8 = newEdge();

    auto f2 = newFace(), f3 = newFace();

    m->position = (v0->position + v1->position)/2;
    m->halfedge() = h0;

    //next, twin, vertex, edge, face
    h0->setNeighbors(h1, h3, m, e0, f0);
    h1->setNeighbors(h6, h1->twin(), v1, e1, f0);
    h2->setNeighbors(h8, h2->twin(), v2, e2, f1);
    h3->setNeighbors(h11, h0, v1, e0, f3);
    h4->setNeighbors(h10, h4->twin(), v0, e3, f2);
    h5->setNeighbors(h3, h5->twin(), v3, e4, f3);
    h6->setNeighbors(h0, h7, v2, e6, f0);
    h7->setNeighbors(h2, h6, m, e6, f1);
    h8->setNeighbors(h7, h9, v0, e7, f1);
    h9->setNeighbors(h4, h8, m, e7, f2);
    h10->setNeighbors(h9, h11, v3, e8, f2);
    h11->setNeighbors(h5, h10, m, e8, f3);

    v0->halfedge() = h8;
    v1->halfedge() = h1;
    v2->halfedge() = h6;
    v3->halfedge() = h10;
    m->halfedge()  = h0;

    e0->halfedge() = h0;
    e6->halfedge() = h6;
    e7->halfedge() = h8;
    e8->halfedge() = h10;

    f0->halfedge() = h0;
    f1->halfedge() = h7;
    f2->halfedge() = h9;
    f3->halfedge() = h11;

    m->isNew = true;
    e0->isNew = false;
    e6->isNew = true;
    e7->isNew = false;
    e8->isNew = true;
    return m;
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.
    HalfedgeIter h;
    VertexIter v;
    for(v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++){
      int n = v->degree();
      double u = n == 3? 3.0/16: 3.0/(8 * n);
      Vector3D sum;
      h = v->halfedge();
      do {
        HalfedgeCIter h_twin = h->twin(); // get the opposite half-edge
        VertexCIter v_n = h_twin->vertex(); // vertex is the 'source' of the half-edge, so
                                          // h->vertex() is v, whereas h_twin->vertex()
                                          // is the neighboring vertex
        sum += v_n->position;
        h = h->twin()->next();               // move to the next outgoing half-edge of the vertex
      } while(h != v->halfedge());          // keep going until we are back where we were
      v->newPosition = (1 - n * u) * v->position + u * sum;
      v->isNew = false;
    }
    printf("phase1:complete\n");
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
    EdgeIter e;
    int num_Edge = 0;
    for(e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++){
      num_Edge++;
      h = e->halfedge();
      auto a = h->vertex();
      auto b = h->twin()->vertex();
      auto c = h->next()->next()->vertex();
      auto d = h->twin()->next()->next()->vertex();
      e->newPosition = ((a->position + b->position) * 3) / 8 + (c->position + d->position) / 8;
      e->isNew = false;
    }
    printf("phase2:complete\n");
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    // e = mesh.edgesBegin();
    // while(e != mesh.edgesEnd()){
    //   EdgeIter nextEdge = e++;
    //   auto pos = e->newPosition;
    //   auto v_new = mesh.splitEdge(e);
    //   v_new->newPosition = pos;
    //   e = nextEdge;
    // }
    // for(e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++){
    //   if(!e->isNew){
    //     auto v_new = mesh.splitEdge(e);
    //     v_new->newPosition = e->newPosition;
    //   }
    //   // auto v_new = mesh.splitEdge(e);
    //   // v_new->newPosition = e->newPosition;
    // }
    e = mesh.edgesBegin();
    for (int i = 0; i < num_Edge; i++, ++e){
      v = mesh.splitEdge(e);
      v->newPosition = e->newPosition;
    }
    printf("phase3:complete\n");
    // 4. Flip any new edge that connects an old and new vertex.
    for(e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++){
      VertexIter s, t;
      h = e->halfedge();
      s = h->vertex();
      t = h->twin()->vertex();
      if(e->isNew)
        if(s->isNew && !t->isNew || !s->isNew && t->isNew)mesh.flipEdge(e);
    }
    // 5. Copy the new vertex positions into final Vertex::position.
    for(v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++){
      v->position = v->newPosition;
    }
  }
}
