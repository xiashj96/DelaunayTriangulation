#include "Delaunay.h"
#include <stdexcept>
#include <list>
#include <algorithm>

namespace dt {
	
	inline double lerp(double start, double end, double t)
	{
		return (1 - t) * start + t * end;
	}

	// 3-rd order determimant
	double determinant(
		double x1, double x2, double x3,
		double y1, double y2, double y3,
		double z1, double z2, double z3)
	{
		return x1 * (y2 * z3 - z2 * y3) + x2 * (y3 * z1 - y1 * z3) + x3 * (y1 * z2 - y2 * z1);
	}

	inline double Area2(const Vertex& p, const Vertex& q, const Vertex& r)
	{
		return determinant(
			1.0, 1.0, 1.0,
			p.x, q.x, r.x,
			p.y, q.y, r.y
		);
	}

	double DistanceSquared(const Vertex& p, const Vertex& q)
	{
		double dx = p.x - q.x;
		double dy = p.y - q.y;
		return dx * dx + dy * dy;
	}

	inline double Distance(const Vertex& p, const Vertex& q)
	{
		return sqrt(DistanceSquared(p, q));
	}


	inline bool ToLeft(const Vertex& p, const Vertex& q, const Vertex& r)
	{
		return Area2(p, q, r) >= 0;
	}

	inline bool Sameline(const Vertex& p, const Vertex& q, const Vertex&   r)
	{
		return Area2(p, q, r) == 0;
	}

	inline bool InTriangle(const Vertex& p, const Vertex& a, const Vertex& b, const Vertex& c)
	{
		return (ToLeft(a, b, p) && ToLeft(b, c, p) && ToLeft(c, a, p));
	}

	// determine if point p is inside the circumcircle of triangle abc
	inline bool InCircle(const Vertex& p, const Vertex& a, const Vertex& b, const Vertex& c)
	{
		double a2 = a.x * a.x + a.y * a.y;
		double b2 = b.x * b.x + b.y * b.y;
		double c2 = c.x * c.x + c.y * c.y;
		double p2 = p.x * p.x + p.y * p.y;

		double d =
			-determinant(
				b.x, c.x, p.x,
				b.y, c.y, p.y,
				b2, c2, p2)
			+ determinant(
				a.x, c.x, p.x,
				a.y, c.y, p.y,
				a2, c2, p2)
			- determinant(
				a.x, b.x, p.x,
				a.y, b.y, p.y,
				a2, b2, p2)
			+ determinant(
				a.x, b.x, c.x,
				a.y, b.y, c.y,
				a2, b2, c2);
		/*double d =
			- determinant(
				b.x, b.y, b2,
				c.x, c.y, c2,
				p.x, p.y, p2)
			+ determinant(
				a.x, a.y, a2,
				c.x, c.y, c2,
				p.x, p.y, p2)
			- determinant(
				a.x, a.y, a2,
				b.x, b.y, b2,
				p.x, p.y, p2)
			+ determinant(
				a.x, a.y, a2,
				b.x, b.y, b2,
				c.x, c.y, c2);*/
		return d > 0.0;
	}

	bool Edge::contains(const Vertex& p) const
	{
		const Vertex& a = *orig;
		const Vertex& b = *(succ->orig);
		double xmax = (a.x > b.x) ? a.x : b.x;
		double xmin = (a.x < b.x) ? a.x : b.x;
		double ymax = (a.y > b.y) ? a.y : b.y;
		double ymin = (a.y < b.y) ? a.y : b.y;
		if (p.x <= xmax && p.x >= xmin && p.y <= ymax && p.y >= ymin)
		{
			if (Sameline(p, a, b))
				return true;
		}
		return false;
	}

	bool Edge::intersect(double height, Vertex& p) const
	{
		const Vertex& p1 = *(this->orig);
		const Vertex& p2 = *(this->succ->orig);
		if ((p1.h > height && p2.h > height) || (p1.h < height && p2.h < height))
		{
			return false;
		}
		
		if (p1.h == height) {
			p = p1;
			return true;
		}

		if (p2.h == height) {
			p = p2;
			return true;
		}

		double t = (height - p1.h) / (p2.h - p1.h);
		p.x = p1.x + t * (p2.x - p1.x);
		p.y = p1.y + t * (p2.y - p1.y);
		p.h = height;
		return true;
	}

	// test if a point is in the triangle (including edges)
	bool Triangle::contains(const Vertex& p) const
	{
		return InTriangle(p, *(inc->orig), *(inc->succ->orig), *(inc->pred->orig));
	}

	bool Triangle::onEdge(const Vertex& p) const
	{
		Edge* ab = inc;
		Edge* bc = inc->succ;
		Edge* ca = inc->pred;
		return ab->contains(p) || bc->contains(p) || ca->contains(p);
	}

	bool Triangle::inCircumcircle(const Vertex& p) const
	{
		return InCircle(p, *(inc->orig), *(inc->succ->orig), *(inc->pred->orig));
	}

	bool Triangle::hasVertex(Vertex* p) const
	{
		if (inc->orig == p || inc->succ->orig == p || inc->pred->orig == p)
		{
			return true;
		}
		else {
			return false;
		}
	}

	void Triangle::interpolate(Vertex& p) const
	{
		if (!contains(p))
		{
			throw std::invalid_argument("query point is not inside the triangle");
		}

		Vertex& a = *(this->inc->orig);
		Vertex& b = *(this->inc->succ->orig);
		Vertex& c = *(this->inc->pred->orig);
		double total = Area2(a, b, c);
		double ta = Area2(b, c, p) / total;
		double tb = Area2(c, a, p) / total;
		double tc = Area2(a, b, p) / total;
		p.h = ta * a.h + tb * b.h + tc * c.h;
		return;
	}

	// rule out wired cases: vertex on plane, edge on plane and triangle on plane
	bool Triangle::intersect(double height) const
	{
		Vertex& a = *(this->inc->orig);
		Vertex& b = *(this->inc->succ->orig);
		Vertex& c = *(this->inc->pred->orig);
		return !((a.h > height && b.h > height && c.h > height) || (a.h <= height && b.h <= height && c.h <= height));
	}

	Delaunay::Delaunay(const std::vector<Vertex>& vertices, double deleteLongTriangleAtBoundary)
	{
		if (vertices.size() < 3)
		{
			throw std::invalid_argument("vertex list is too small");
		}
		
		maxX = minX = vertices[0].x;
		maxY = minY = vertices[0].y;
		maxH = minH = vertices[0].h;

		for (const Vertex& v : vertices)
		{
			if (v.x > maxX) maxX = v.x;
			if (v.x < minX) minX = v.x;
			if (v.y > maxY) maxY = v.y;
			if (v.y < minY) minY = v.y;
			if (v.h > maxH) maxH = v.h;
			if (v.h < minH) minH = v.h;
		}

		// normalize and store the vertices locally
		_vertices = vertices;
		for (Vertex& v : _vertices)
		{
			v = world2screen(v);
		}
		triangulate(deleteLongTriangleAtBoundary);
	}

	Delaunay::~Delaunay()
	{
		for (const auto& e : _edges)
		{
			delete e;
		}
		for (const auto& t : _triangles)
		{
			delete t;
		}
	}

	void Delaunay::triangulate(double deleteLongTriangleAtBoundary)
	{
		// initilization: add super triangle
		double size = lerp(3.0, 1.1, deleteLongTriangleAtBoundary);
		Vertex a(-size, -size);
		Vertex b(size, -size);
		Vertex c(size, size);
		Vertex d(-size, size);

		Edge* ab = new Edge;
		Edge* bc = new Edge;
		Edge* ca = new Edge;
		Edge* ac = new Edge;
		Edge* cd = new Edge;
		Edge* da = new Edge;
		Triangle* abc = new Triangle;
		Triangle* acd = new Triangle;
		
		ab->orig = &a;
		ab->pred = ca;
		ab->succ = bc;
		ab->twin = nullptr;
		ab->inc = abc;

		bc->orig = &b;
		bc->pred = ab;
		bc->succ = ca;
		bc->twin = nullptr;
		bc->inc = abc;

		ca->orig = &c;
		ca->pred = bc;
		ca->succ = ab;
		ca->twin = ac;
		ca->inc = abc;

		ac->orig = &a;
		ac->pred = da;
		ac->succ = cd;
		ac->twin = ca;
		ac->inc = acd;

		cd->orig = &c;
		cd->pred = ac;
		cd->succ = da;
		cd->twin = nullptr;
		cd->inc = acd;

		da->orig = &d;
		da->pred = cd;
		da->succ = ac;
		da->twin = nullptr;
		da->inc = acd;

		abc->inc = ab;
		acd->inc = ac;

		_edges.push_back(ab);
		_edges.push_back(bc);
		_edges.push_back(ca);
		_edges.push_back(ac);
		_edges.push_back(cd);
		_edges.push_back(da);
		_triangles.push_back(abc);
		_triangles.push_back(acd);

		int index = 0;

		// add new vertices in random order
		for (auto& v : _vertices)
		{
			// point location: determine which triangle v is in
			index++;
			Triangle* t = nullptr;
			for (const auto& triangle : _triangles)
			{
				if (triangle->contains(v))
				{
					t = triangle;
					break;
				}
			}
			if (!t)
			{
				throw std::logic_error("new vertex is not contained in an existing triangle");
				// this should never happen, since we initialize _triangles with a super triangle
			}
			else {
				// add 2 new triangles
				if (!t->onEdge(v)) // general case: p is strictly inside a triangle
				{
					// 1. construct 3 new edges (6 half-edges) and 2 new triangles
					Vertex* a = t->inc->orig;
					Vertex* b = t->inc->succ->orig;
					Vertex* c = t->inc->pred->orig;

					Edge* ab = t->inc;
					Edge* bc = t->inc->succ;
					Edge* ca = t->inc->pred;
					Edge* va = new Edge;
					Edge* vb = new Edge;
					Edge* vc = new Edge;
					Edge* av = new Edge;
					Edge* bv = new Edge;
					Edge* cv = new Edge;

					Triangle* abv = t;
					Triangle* bcv = new Triangle;
					Triangle* cav = new Triangle;

					ab->pred = va;
					ab->succ = bv;
					ab->inc = abv;

					bc->pred = vb;
					bc->succ = cv;
					bc->inc = bcv;

					ca->pred = vc;
					ca->succ = av;
					ca->inc = cav;

					va->orig = &v;
					va->pred = bv;
					va->succ = ab;
					va->twin = av;
					va->inc = abv;

					av->orig = a;
					av->pred = ca;
					av->succ = vc;
					av->twin = va;
					av->inc = cav;

					vb->orig = &v;
					vb->pred = cv;
					vb->succ = bc;
					vb->twin = bv;
					vb->inc = bcv;

					bv->orig = b;
					bv->pred = ab;
					bv->succ = va;
					bv->twin = vb;
					bv->inc = abv;

					vc->orig = &v;
					vc->pred = av;
					vc->succ = ca;
					vc->twin = cv;
					vc->inc = cav;

					cv->orig = c;
					cv->pred = bc;
					cv->succ = vb;
					cv->twin = vc;
					cv->inc = bcv;

					abv->inc = ab; // unnecessary
					bcv->inc = bc;
					cav->inc = ca;

					//v.inc = va;
					//a->inc = ab;
					//b->inc = bc;
					//c->inc = ca;

					// legalize all edges by edge flipping
					legalizeEdge(&v, ab);
					legalizeEdge(&v, bc);
					legalizeEdge(&v, ca);

					_edges.push_back(av);
					_edges.push_back(bv);
					_edges.push_back(cv);
					_edges.push_back(va);
					_edges.push_back(vb);
					_edges.push_back(vc);

					_triangles.push_back(bcv);
					_triangles.push_back(cav);

				}
				else {	// special case: p is on an edge
					// figure out which edge it is on and make it t->inc
					while (!t->inc->contains(v))
					{
						t->inc = t->inc->succ;
					}

					if (t->inc->isExtreme()) // v should not be on an extreme edge since we are working inside a super triangle
					{
						throw std::logic_error("new point on extreme edge");
					}

					// add 2 new triangles and 6 new half edges
					Vertex* a = t->inc->orig;
					Vertex* b = t->inc->succ->orig;
					Vertex* c = t->inc->pred->orig;
					Vertex* d = t->inc->twin->pred->orig;

					Edge* bc = t->inc->succ;
					Edge* ca = t->inc->pred;
					Edge* ad = t->inc->twin->succ;
					Edge* db = t->inc->twin->pred;
					Edge* av = t->inc;
					Edge* bv = t->inc->twin;
					Edge* cv = new Edge;
					Edge* dv = new Edge;
					Edge* va = new Edge;
					Edge* vb = new Edge;
					Edge* vc = new Edge;
					Edge* vd = new Edge;

					Triangle* avc = t;
					Triangle* bvd = t->inc->twin->inc;
					Triangle* cvb = new Triangle;
					Triangle* dva = new Triangle;

					// update edges
					bc->succ = cv;
					bc->pred = vb;
					bc->inc = cvb;

					ca->succ = av;
					ca->pred = vc;
					ca->inc = avc;

					ad->succ = dv;
					ad->pred = va;
					ad->inc = dva;

					db->succ = bv;
					db->pred = vd;
					db->inc = bvd;

					av->orig = a;	// unnecessary
					av->succ = vc;
					av->pred = ca;	// unnecessary
					av->twin = va;
					av->inc = avc;	// unnecessary

					bv->orig = b;
					bv->succ = vd;
					bv->pred = db;
					bv->twin = vb;
					bv->inc = bvd;

					cv->orig = c;
					cv->succ = vb;
					cv->pred = bc;
					cv->twin = vc;
					cv->inc = cvb;

					dv->orig = d;
					dv->succ = va;
					dv->pred = ad;
					dv->twin = vd;
					dv->inc = dva;

					va->orig = &v;
					va->pred = dv;
					va->succ = ad;
					va->twin = av;
					va->inc = dva;

					vb->orig = &v;
					vb->succ = bc;
					vb->pred = cv;
					vb->twin = bv;
					vb->inc = cvb;

					vc->orig = &v;
					vc->succ = ca;
					vc->pred = av;
					vc->twin = cv;
					vc->inc = avc;

					vd->orig = &v;
					vd->succ = db;
					vd->pred = bv;
					vd->twin = dv;
					vd->inc = bvd;

					// update triangles
					avc->inc = av;
					cvb->inc = cv;
					bvd->inc = bv;
					dva->inc = dv;
					
					// update vertex
					//v.inc = va;
					//a->inc = ad;
					//b->inc = bc;
					//c->inc = ca;
					//d->inc = db;

					// legalize all edges by edge flipping
					legalizeEdge(&v, ad);
					legalizeEdge(&v, db);
					legalizeEdge(&v, bc);
					legalizeEdge(&v, ca);

					_edges.push_back(cv);
					_edges.push_back(dv);
					_edges.push_back(va);
					_edges.push_back(vb);
					_edges.push_back(vc);
					_edges.push_back(vd);

					_triangles.push_back(cvb);
					_triangles.push_back(dva);
				}
			}
		}

		// delete all triangles connected to 3 points of the super triangle
		for (auto& t : _triangles)
		{
			if (t->hasVertex(&a) || t->hasVertex(&b) || t->hasVertex(&c) || t->hasVertex(&d))
			{
				// also make 3 incident edges' incident triangle null
				t->inc->inc = nullptr;
				t->inc->pred->inc = nullptr;
				t->inc->succ->inc = nullptr;
				delete t;
				t = nullptr;
			}
		}
		// delete corresponding edges
		for (auto& e : _edges)
		{
			if (e->inc == nullptr)
			{
				// also make sure their twins edges(if exist) has empty twin
				if (!e->isExtreme())
				{
					e->twin->twin = nullptr;
				}
				delete e;
				e = nullptr;
			}
		}

		_triangles.erase(std::remove(_triangles.begin(), _triangles.end(), nullptr), _triangles.end());
		_edges.erase(std::remove(_edges.begin(), _edges.end(), nullptr), _edges.end());

		return;
	}

	std::vector<unsigned int> Delaunay::getTriangleVertexIndices() const
	{
		std::vector<unsigned int> indices;
		for (const auto& t : _triangles)
		{
			unsigned int a = t->inc->orig - _vertices.data();
			unsigned int b = t->inc->succ->orig - _vertices.data();
			unsigned int c = t->inc->pred->orig - _vertices.data();
			indices.push_back(a);
			indices.push_back(b);
			indices.push_back(c);
		}
		return indices;
	}

	bool Delaunay::interpolate(Vertex& p_world) const
	{
		// point location: in which triangle is p located?
		Vertex p_screen = world2screen(p_world);
		Triangle* t = nullptr;
		for (const auto& triangle : _triangles)
		{
			if (triangle->contains(p_screen))
			{
				t = triangle;
				break;
			}
		}
		if (!t)
		{
			return false;
		}
		else {
			t->interpolate(p_screen);
			p_world.h = screen2world(p_screen).h;
			return true;
		}
	}

	std::vector<Contour> Delaunay::getContours(double height)
	{
		std::vector<Contour> contours;
		if (height < minH || height > maxH)
		{
			return contours;
		}
		height = (height - minH) / (maxH - minH) * 2 - 1; // world2screen

		std::vector<Triangle*> unvisited = _triangles;
		while (!unvisited.empty()) // end condition: all triangles are visited
		{
			// get a random triangle in the unvisited collection
			Triangle* t = unvisited.back();

			// if t intersects height, start tracing a contour line in both directions
			if (!t->intersect(height))
			{
				unvisited.pop_back();
			}
			else {
				std::list<Vertex> c;
				Vertex v;
				Edge* e = t->inc;
				// find the first edge that intersects h=height and also has orig->h > h;
				while (!((e->intersect(height, v)) && (e->orig->h > height)))
				{
					e = e->succ;
				}
				Edge* start_edge = e;

				c.push_back(v);
				// every interation process a triangle along the contour
				// until it reaches the boundary or back to the start point (forming a closed loop)
				// for h = some point exact height, this algorithm will create some duplicate points in the contour, postprocess them to remove the duplcate points
				while (e != nullptr)
				{
					auto cur = std::find(unvisited.begin(), unvisited.end(), e->inc);
					if (cur == unvisited.end()) break;
					if (e->pred->intersect(height, v))
					{
						c.push_back(v);
						e = e->pred->twin;
					}
					else if (e->succ->intersect(height, v))
					{
						c.push_back(v);
						e = e->succ->twin;
					}
					unvisited.erase(cur);
				};

				// the other direction
				e = start_edge->twin;
				while (e != nullptr)
				{
					auto cur = std::find(unvisited.begin(), unvisited.end(), e->inc);
					if (cur == unvisited.end()) break;
					if (e->succ->intersect(height, v))
					{ 
						c.push_front(v);
						e = e->succ->twin;
					}
					else if (e->pred->intersect(height, v))
					{
						c.push_front(v);
						e = e->pred->twin;
					}
					unvisited.erase(cur);
				}
				c.unique(); // remove consecutive duplicate points

				if (c.size() > 1)
				{
					for (auto& v : c)
					{
						v = screen2world(v);
					}
					contours.emplace_back(c.begin(), c.end());
				}
			}
		}
		return contours;
	}

	Vertex Delaunay::world2screen(const Vertex& p) const
	{
		double x = (p.x - minX) / (maxX - minX) * 2 - 1;
		double y = (p.y - minY) / (maxY - minY) * 2 - 1;
		double h = (p.h - minH) / (maxH - minH) * 2 - 1;
		return Vertex(x, y, h);
	}

	Vertex Delaunay::screen2world(const Vertex& p) const
	{
		double x = (p.x + 1) / 2 * (maxX - minX) + minX;
		double y = (p.y + 1) / 2 * (maxY - minY) + minY;
		double h = (p.h + 1) / 2 * (maxH - minH) + minH;
		return Vertex(x, y, h);
	}

	void Delaunay::legalizeEdge(Vertex* p, Edge* e)
	{
		if (e->isExtreme())
		{
			return;
		}

		// assume pab is a newly created CCW triangle inside frontier
		Vertex* a = e->orig;
		Vertex* b = e->twin->orig;
		Vertex* q = e->twin->pred->orig;
		Edge* pa = e->pred;
		Edge* bp = e->succ;
		Edge* qb = e->twin->pred;
		Edge* aq = e->twin->succ;
		Triangle* t1 = e->inc;
		Triangle* t2 = e->twin->inc;

		//if (q->h == 0) return;

		if (t1->inCircumcircle(*q))
		{
			// perform edge flip
			// t1: abp->pqb, t2: baq->qpa
			// update edges
			Edge* pq = e;
			Edge* qp = e->twin;

			pa->orig = p;
			pa->succ = aq;
			pa->pred = qp;
			pa->inc = t2;

			bp->orig = b;
			bp->succ = pq;
			bp->pred = qb;
			bp->inc = t1;

			qb->orig = q;
			qb->succ = bp;
			qb->pred = pq;
			qb->inc = t1;

			aq->orig = a;
			aq->succ = qp;
			aq->pred = pa;
			aq->inc = t2;

			pq->orig = p;
			pq->succ = qb;
			pq->pred = bp;
			pq->inc = t1;

			qp->orig = q;
			qp->succ = pa;
			qp->pred = aq;
			qp->inc = t2;

			// update triangles
			t1->inc = pq;
			t2->inc = qp;

			// recursive call
			legalizeEdge(p, qb);
			legalizeEdge(p, aq);
		}
		return;
	}


}