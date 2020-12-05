#pragma once
#include <vector>

namespace dt {
	
	/*
		given a 2d point set P, solve the following problems:
		1. construct a Delaunay triangulation (DT)
		2. given ant location inside the convex hull of P, calculate linear interpolated height value based on DT
		3. generate piecewise linear contour line from DT
		4. smoothed height interpolation based on DT
		5. smoothed contour lines

		the complexity of the algorithm should not exceed O(nlogn)

		TODO:
		- color mapping fragment shader
		- contour line tracing algorithm
		- contour drawing by GL_LINE_STRIP and GL_LINE_LOOP
		- font rendering, display x,y and h values in left bottom
	*/

	struct Vertex;
	struct Edge;
	struct Triangle;

	// to allow efficient tracing of contour lines among the triangle network, we use DCEL to store the topological structure
	struct Vertex {
		Vertex(double vx = 0.0, double vy = 0.0, double vh = 0.0) : x(vx), y(vy), h(vh) {}
		double x;
		double y;
		double h; // height
		bool operator==(const Vertex& other) const { return x == other.x && y == other.y; }
		bool operator!=(const Vertex& other) const { return !(*this == other); }
		bool operator<(const Vertex& other) const { // lexicographical order
			if (y < other.y)
				return true;
			else if (y == other.y)
				return x < other.x;
			else
				return false;
		}
			
		//Edge* inc = nullptr; // "first" outgoing half edge
	};

	struct Edge {
		Vertex* orig = nullptr;
		Edge* pred = nullptr;
		Edge* succ = nullptr;
		Triangle* inc = nullptr; // left triangle
		Edge* twin = nullptr;
		bool isExtreme() const
		{
			return twin == nullptr;
		}
		bool contains(const Vertex& p) const;
		bool intersect(double height, Vertex& p) const;
	};

	struct Triangle {
		Edge* inc = nullptr; // "first" incident edge
		bool contains(const Vertex& p) const;
		bool onEdge(const Vertex& p) const;
		bool inCircumcircle(const Vertex& p) const;
		bool hasVertex(Vertex* p) const;
		void interpolate(Vertex& p) const;
		bool intersect(double height) const;
	};

	using Contour = std::vector<Vertex>;

	class Delaunay
	{
	public:
		Delaunay(const std::vector<Vertex>& vertices);
		Delaunay(const Delaunay&) = delete;
		Delaunay(Delaunay&&) = delete;
		Delaunay& operator=(const Delaunay&) = delete;
		Delaunay& operator=(Delaunay&&) = delete;
		~Delaunay();

		unsigned int getNumOfTriangles() const { return _triangles.size(); }
		const std::vector<Vertex>& getVertices() const { return _vertices; }
		std::vector<unsigned int> getTriangleVertexIndices() const;
		// return false if p is outside the convex hull of point set
		bool interpolate(Vertex& p) const;
		std::vector<Contour> getContours(double height, double minRelativeLength = 0.1);

		double getMaxX() const { return maxX; }
		double getMinX() const { return minX; }
		double getMaxY() const { return maxY; }
		double getMinY() const { return minY; }
		double getMaxH() const { return maxH; }
		double getMinH() const { return minH; }
	private:
		void triangulate();
		void legalizeEdge(Vertex* p, Edge* e);

		double minX, maxX, minY, maxY, minH, maxH;
		// note: triangles and edges are stored as pointers, need to dispose memory in destructor
		std::vector<Triangle*> _triangles;
		std::vector<Edge*> _edges;
		std::vector<Vertex> _vertices;
	};



}

