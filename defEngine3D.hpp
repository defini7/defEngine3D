#pragma once

/*
* TODO:
* - Lighting
* - Mipmapping
* - Texture sampling
* - Perspective correction !!!
*/

#include "defGameEngine.h"

#define DEF_MATH_3D
#include "defMath3D.hpp"

#pragma warning(disable : 4996)

#include <fstream>
#include <strstream>

namespace def::gfx3d
{
	template <class T>
	struct Camera3
	{
		Camera3() = default;

		math::Vector3<T> pos;
		math::Vector3<T> dir;
		math::Vector3<T> rot;

		void rotate_by(const math::Vector3<T>& v);
		void move_by(const math::Vector3<T>& v);
	};

	typedef Camera3<int> Camera3i;
	typedef Camera3<float> Camera3f;
	typedef Camera3<double> Camera3d;

	template <class T>
	struct Triangle
	{
		math::Vector3<T> p[3];
		math::Tex3<T> t[3];
		def::Pixel col[3];

		uint32_t clip_against_plane(math::Vector3<T> planeP, math::Vector3<T> planeN, Triangle<T>& outTri1, Triangle<T>& outTri2);
	};

	typedef Triangle<int> Triangle3i;
	typedef Triangle<float> Triangle3f;
	typedef Triangle<double> Triangle3d;

	struct Mesh
	{
		std::vector<Triangle3f> tris;

		bool load_obj(const std::string& objFilename, bool textured = false);
		
		static Mesh MakeCube();
	};

	class Engine
	{
	public:
		Engine(
			const def::vi2d& viewPort,
			float fovRad,
			float aspectRatio,
			float near = 0.001f,
			float far = 1000.0f
		);

		~Engine();

	private:
		math::Matrix m_Proj;
		math::Matrix m_WorldMat;

		float* m_DepthBuffer = nullptr;

		def::vi2d m_ViewPort;

		Mesh m_Model;

		def::Graphic* m_DrawTarget = nullptr;
		const def::Sprite* m_TargetTexture = nullptr;
		def::Pixel m_TargetCol[3] = { def::WHITE, def::WHITE, def::WHITE };

		std::list<Triangle3f> m_TrianglesToRaster;

		Triangle3f* m_FormingTriangle;

		size_t m_Vertex;
		size_t m_TexCoord;

		Sprite::SampleMethod m_TexSampleMethod;
		Sprite::WrapMethod m_TexWrapMethod;

	public:
		def::Graphic* DrawTarget();
		math::Matrix& WorldMatrix();

		Sprite::SampleMethod& SampleMethod();
		Sprite::WrapMethod& WrapMethod();

		void Begin();
		void End();

		void Texture(const def::Sprite* texture);
		void Color(const def::Pixel& p1, const def::Pixel& p2, const def::Pixel& p3);
		void Vertex(const math::Vector3f& v);
		void TexCoord(const math::Tex3f& t);
		void SetModel(const Mesh& model);

		void Clear(const def::Pixel& col, bool clearDrawTarget = true, bool clearDepthBuffer = true);
		void Update(Camera3f& camera, const math::Vector3f& up = { 0.0f, 1.0f, 0.0f }, const math::Vector3f& target = { 0.0f, 0.0f, 1.0f });
		void Draw();

		void DrawTexturedTriangle(
			def::vi2d p1, math::Tex3f t1,
			def::vi2d p2, math::Tex3f t2,
			def::vi2d p3, math::Tex3f t3,
			const def::Pixel& col1,
			const def::Pixel& col2,
			const def::Pixel& col3,
			const def::Sprite* spr
		);

	};

#ifdef DEF_ENGINE_3D
#undef DEF_ENGINE_3D

	template <class T> void Camera3<T>::rotate_by(const math::Vector3<T>& v) { rot += v; }
	template <class T> void Camera3<T>::move_by(const math::Vector3<T>& v) { pos += v; }

	template <class T>
	uint32_t Triangle<T>::clip_against_plane(math::Vector3<T> planeP, math::Vector3<T> planeN, Triangle<T>& outTri1, Triangle<T>& outTri2)
	{
		planeN = planeN.norm();

		auto dist = [&](math::Vector3<T>& p)
		{
			return planeN.dot(p) - planeN.dot(planeP);
		};

		uint32_t insidePointsCount = 0;
		uint32_t outsidePointsCount = 0;

		math::Vector3<T>* insidePoints[3];
		math::Vector3<T>* outsidePoints[3];

		uint32_t insideTexturesCount = 0;
		uint32_t outsideTexturesCount = 0;

		math::Tex3<T>* insideTextures[3];
		math::Tex3<T>* outsideTextures[3];

		float d0 = dist(p[0]);
		float d1 = dist(p[1]);
		float d2 = dist(p[2]);

		if (d0 >= 0)
		{
			insidePoints[insidePointsCount++] = &p[0];
			insideTextures[insideTexturesCount++] = &t[0];
		}
		else
		{
			outsidePoints[outsidePointsCount++] = &p[0];
			outsideTextures[outsideTexturesCount++] = &t[0];
		}

		if (d1 >= 0)
		{
			insidePoints[insidePointsCount++] = &p[1];
			insideTextures[insideTexturesCount++] = &t[1];
		}
		else
		{
			outsidePoints[outsidePointsCount++] = &p[1];
			outsideTextures[outsideTexturesCount++] = &t[1];
		}

		if (d2 >= 0)
		{
			insidePoints[insidePointsCount++] = &p[2];
			insideTextures[insideTexturesCount++] = &t[2];
		}
		else
		{
			outsidePoints[outsidePointsCount++] = &p[2];
			outsideTextures[outsideTexturesCount++] = &t[2];
		}

		if (insidePointsCount == 0) return 0;

		if (insidePointsCount == 3)
		{
			outTri1 = *this;
			return 1;
		}

		if (insidePointsCount == 1 && outsidePointsCount == 2)
		{
			for (size_t i = 0; i < 3; i++) outTri1.col[i] = col[i];

			outTri1.p[0] = *insidePoints[0];
			outTri1.t[0] = *insideTextures[0];

			float t;

			outTri1.p[1] = insidePoints[0]->intersect_plane(*outsidePoints[0], planeP, planeN, t);
			outTri1.t[1] = t * (*outsideTextures[0] - *insideTextures[0]) + *insideTextures[0];

			outTri1.p[2] = insidePoints[0]->intersect_plane(*outsidePoints[1], planeP, planeN, t);
			outTri1.t[2] = t * (*outsideTextures[1] - *insideTextures[0]) + *insideTextures[0];

			return 1;
		}

		if (insidePointsCount == 2 && outsidePointsCount == 1)
		{
			for (size_t i = 0; i < 3; i++)
			{
				outTri1.col[i] = col[i];
				outTri2.col[i] = col[i];
			}

			outTri1.p[0] = *insidePoints[0];
			outTri1.p[1] = *insidePoints[1];

			outTri1.t[0] = *insideTextures[0];
			outTri1.t[1] = *insideTextures[1];

			float t;

			outTri1.p[2] = insidePoints[0]->intersect_plane(*outsidePoints[0], planeP, planeN, t);
			outTri1.t[2] = t * (*outsideTextures[0] - *insideTextures[0]) + *insideTextures[0];

			outTri2.p[0] = *insidePoints[1];
			outTri2.t[0] = *insideTextures[1];

			outTri2.p[1] = outTri1.p[2];
			outTri2.t[1] = outTri1.t[2];

			outTri2.p[2] = insidePoints[1]->intersect_plane(*outsidePoints[0], planeP, planeN, t);
			outTri2.t[2] = t * (*outsideTextures[0] - *insideTextures[1]) + *insideTextures[1];

			return 2;
		}
	}


	bool Mesh::load_obj(const std::string& objFilename, bool textured)
	{
		std::ifstream f(objFilename);
		if (!f.is_open()) return false;

		std::vector<math::Vector3f> verts;
		std::vector<math::Tex3f> texs;

		while (!f.eof())
		{
			char line[128];
			f.getline(line, 128);

			std::strstream s;
			s << line;

			char junk;

			if (line[0] == 'v')
			{
				if (line[1] == 't')
				{
					math::Tex3f t;
					s >> junk >> junk >> t.u >> t.v >> t.w;
					t.v = 1.0f - t.v;
					texs.push_back(t);
				}
				else
				{
					math::Vector3f p;
					s >> junk >> p.x >> p.y >> p.z;
					verts.push_back(p);
				}
			}

			if (textured)
			{
				if (line[0] == 'f')
				{
					s >> junk;

					std::string tokens[6];
					int tokenCount = -1;
					while (!s.eof())
					{
						char c = s.get();
						if (c == ' ' || c == '/')
						{
							if (tokenCount > 0)
							{
								if (!tokens[tokenCount].empty())
									tokenCount++;
							}
							else
								tokenCount++;
						}
						else
							tokens[tokenCount].append(1, c);
					}

					tokens[tokenCount].pop_back();

					tris.push_back({
						{
							verts[stoi(tokens[0]) - 1],
							verts[stoi(tokens[2]) - 1],
							verts[stoi(tokens[4]) - 1]
						},
						{
							texs[stoi(tokens[1]) - 1],
							texs[stoi(tokens[3]) - 1],
							texs[stoi(tokens[5]) - 1],
						}
					});
				}
			}
			else
			{
				if (line[0] == 'f')
				{
					int f[3];
					s >> junk >> f[0] >> f[1] >> f[2];
					tris.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
				}
			}
		}

		return true;
	}

	Mesh Mesh::MakeCube()
	{
		Mesh m;

		m.tris =
		{
#ifdef DEMO_TEX
			// SOUTH 
			def::gfx3d::Triangle3f{ def::math::Vector3{ 0.0f, 0.0f, 0.0f, 1.0f },    def::math::Vector3{ 0.0f, 1.0f, 0.0f, 1.0f },    def::math::Vector3{ 1.0f, 1.0f, 0.0f, 1.0f },		def::math::Tex3{ -1.0f, 2.0f, 1.0f },		def::math::Tex3{ -1.0f, -1.0f, 1.0f },		def::math::Tex3{ 2.0f, -1.0f, 1.0f } },
			def::gfx3d::Triangle3f{ def::math::Vector3{ 0.0f, 0.0f, 0.0f, 1.0f },    def::math::Vector3{ 1.0f, 1.0f, 0.0f, 1.0f },    def::math::Vector3{ 1.0f, 0.0f, 0.0f, 1.0f },		def::math::Tex3{ -1.0f, 2.0f, 1.0f },		def::math::Tex3{ 2.0f, -1.0f, 1.0f },		def::math::Tex3{ 2.0f, 2.0f, 1.0f } },

			// EAST           																			   
			def::gfx3d::Triangle3f{ def::math::Vector3{ 1.0f, 0.0f, 0.0f, 1.0f },    def::math::Vector3{ 1.0f, 1.0f, 0.0f, 1.0f },    def::math::Vector3{ 1.0f, 1.0f, 1.0f, 1.0f },		def::math::Tex3{ -1.0f, 2.0f, 1.0f },		def::math::Tex3{ -1.0f, -1.0f, 1.0f },		def::math::Tex3{ 2.0f, -1.0f, 1.0f } },
			def::gfx3d::Triangle3f{ def::math::Vector3{ 1.0f, 0.0f, 0.0f, 1.0f },    def::math::Vector3{ 1.0f, 1.0f, 1.0f, 1.0f },    def::math::Vector3{ 1.0f, 0.0f, 1.0f, 1.0f },		def::math::Tex3{ -1.0f, 2.0f, 1.0f },		def::math::Tex3{ 2.0f, -1.0f, 1.0f },		def::math::Tex3{ 2.0f, 2.0f, 1.0f } },

			// NORTH           																			   
			def::gfx3d::Triangle3f{ def::math::Vector3{ 1.0f, 0.0f, 1.0f, 1.0f },    def::math::Vector3{ 1.0f, 1.0f, 1.0f, 1.0f },    def::math::Vector3{ 0.0f, 1.0f, 1.0f, 1.0f },		def::math::Tex3{ -1.0f, 2.0f, 1.0f },		def::math::Tex3{ -1.0f, -1.0f, 1.0f},		def::math::Tex3{ 2.0f, -1.0f, 1.0f } },
			def::gfx3d::Triangle3f{ def::math::Vector3{ 1.0f, 0.0f, 1.0f, 1.0f },    def::math::Vector3{ 0.0f, 1.0f, 1.0f, 1.0f },    def::math::Vector3{ 0.0f, 0.0f, 1.0f, 1.0f },		def::math::Tex3{ -1.0f, 2.0f, 1.0f },		def::math::Tex3{ 2.0f, -1.0f, 1.0f},		def::math::Tex3{ 2.0f, 2.0f, 1.0f } },

			// WEST            																			   
			def::gfx3d::Triangle3f{ def::math::Vector3{ 0.0f, 0.0f, 1.0f, 1.0f },    def::math::Vector3{ 0.0f, 1.0f, 1.0f, 1.0f },    def::math::Vector3{ 0.0f, 1.0f, 0.0f, 1.0f },		def::math::Tex3{ -1.0f, 2.0f, 1.0f },		def::math::Tex3{ -1.0f, -1.0f, 1.0f },		def::math::Tex3{ 2.0f, -1.0f, 1.0f } },
			def::gfx3d::Triangle3f{ def::math::Vector3{ 0.0f, 0.0f, 1.0f, 1.0f },    def::math::Vector3{ 0.0f, 1.0f, 0.0f, 1.0f },    def::math::Vector3{ 0.0f, 0.0f, 0.0f, 1.0f },		def::math::Tex3{ -1.0f, 2.0f, 1.0f },		def::math::Tex3{ 2.0f, -1.0f, 1.0f },		def::math::Tex3{ 2.0f, 2.0f, 1.0f } },

			// TOP             																			   
			def::gfx3d::Triangle3f{ def::math::Vector3{ 0.0f, 1.0f, 0.0f, 1.0f },    def::math::Vector3{ 0.0f, 1.0f, 1.0f, 1.0f },    def::math::Vector3{ 1.0f, 1.0f, 1.0f, 1.0f },		def::math::Tex3{ -1.0f, 2.0f, 1.0f },		def::math::Tex3{ -1.0f, -1.0f, 1.0f },		def::math::Tex3{ 2.0f, -1.0f, 1.0f } },
			def::gfx3d::Triangle3f{ def::math::Vector3{ 0.0f, 1.0f, 0.0f, 1.0f },    def::math::Vector3{ 1.0f, 1.0f, 1.0f, 1.0f },    def::math::Vector3{ 1.0f, 1.0f, 0.0f, 1.0f },		def::math::Tex3{ -1.0f, 2.0f, 1.0f },		def::math::Tex3{ 2.0f, -1.0f, 1.0f },		def::math::Tex3{ 2.0f, 2.0f, 1.0f } },

			// BOTTOM          																			  
			def::gfx3d::Triangle3f{ def::math::Vector3{ 1.0f, 0.0f, 1.0f, 1.0f },    def::math::Vector3{ 0.0f, 0.0f, 1.0f, 1.0f },    def::math::Vector3{ 0.0f, 0.0f, 0.0f, 1.0f },		def::math::Tex3{ -1.0f, 2.0f, 1.0f },		def::math::Tex3{ -1.0f, -1.0f, 1.0f },		def::math::Tex3{ 2.0f, -1.0f, 1.0f } },
			def::gfx3d::Triangle3f{ def::math::Vector3{ 1.0f, 0.0f, 1.0f, 1.0f },    def::math::Vector3{ 0.0f, 0.0f, 0.0f, 1.0f },    def::math::Vector3{ 1.0f, 0.0f, 0.0f, 1.0f },		def::math::Tex3{ -1.0f, 2.0f, 1.0f },		def::math::Tex3{ 2.0f, -1.0f, 1.0f },		def::math::Tex3{ 2.0f, 2.0f, 1.0f } },
#else
			// SOUTH
			def::gfx3d::Triangle3f{ def::math::Vector3{ 0.0f, 0.0f, 0.0f, 1.0f },    def::math::Vector3{ 0.0f, 1.0f, 0.0f, 1.0f },    def::math::Vector3{ 1.0f, 1.0f, 0.0f, 1.0f },		def::math::Tex3{ 0.0f, 1.0f, 1.0f },		def::math::Tex3{ 0.0f, 0.0f, 1.0f },		def::math::Tex3{ 1.0f, 0.0f, 1.0f } },
			def::gfx3d::Triangle3f{ def::math::Vector3{ 0.0f, 0.0f, 0.0f, 1.0f },    def::math::Vector3{ 1.0f, 1.0f, 0.0f, 1.0f },    def::math::Vector3{ 1.0f, 0.0f, 0.0f, 1.0f },		def::math::Tex3{ 0.0f, 1.0f, 1.0f },		def::math::Tex3{ 1.0f, 0.0f, 1.0f },		def::math::Tex3{ 1.0f, 1.0f, 1.0f } },

			// EAST           																			   
			def::gfx3d::Triangle3f{ def::math::Vector3{ 1.0f, 0.0f, 0.0f, 1.0f },    def::math::Vector3{ 1.0f, 1.0f, 0.0f, 1.0f },    def::math::Vector3{ 1.0f, 1.0f, 1.0f, 1.0f },		def::math::Tex3{ 0.0f, 1.0f, 1.0f },		def::math::Tex3{ 0.0f, 0.0f, 1.0f },		def::math::Tex3{ 1.0f, 0.0f, 1.0f } },
			def::gfx3d::Triangle3f{ def::math::Vector3{ 1.0f, 0.0f, 0.0f, 1.0f },    def::math::Vector3{ 1.0f, 1.0f, 1.0f, 1.0f },    def::math::Vector3{ 1.0f, 0.0f, 1.0f, 1.0f },		def::math::Tex3{ 0.0f, 1.0f, 1.0f },		def::math::Tex3{ 1.0f, 0.0f, 1.0f },		def::math::Tex3{ 1.0f, 1.0f, 1.0f } },

			// NORTH           																			   
			def::gfx3d::Triangle3f{ def::math::Vector3{ 1.0f, 0.0f, 1.0f, 1.0f },    def::math::Vector3{ 1.0f, 1.0f, 1.0f, 1.0f },    def::math::Vector3{ 0.0f, 1.0f, 1.0f, 1.0f },		def::math::Tex3{ 0.0f, 1.0f, 1.0f },		def::math::Tex3{ 0.0f, 0.0f, 1.0f},			def::math::Tex3{ 1.0f, 0.0f, 1.0f } },
			def::gfx3d::Triangle3f{ def::math::Vector3{ 1.0f, 0.0f, 1.0f, 1.0f },    def::math::Vector3{ 0.0f, 1.0f, 1.0f, 1.0f },    def::math::Vector3{ 0.0f, 0.0f, 1.0f, 1.0f },		def::math::Tex3{ 0.0f, 1.0f, 1.0f },		def::math::Tex3{ 1.0f, 0.0f, 1.0f},			def::math::Tex3{ 1.0f, 1.0f, 1.0f } },

			// WEST            																			   
			def::gfx3d::Triangle3f{ def::math::Vector3{ 0.0f, 0.0f, 1.0f, 1.0f },    def::math::Vector3{ 0.0f, 1.0f, 1.0f, 1.0f },    def::math::Vector3{ 0.0f, 1.0f, 0.0f, 1.0f },		def::math::Tex3{ 0.0f, 1.0f, 1.0f },		def::math::Tex3{ 0.0f, 0.0f, 1.0f },		def::math::Tex3{ 1.0f, 0.0f, 1.0f } },
			def::gfx3d::Triangle3f{ def::math::Vector3{ 0.0f, 0.0f, 1.0f, 1.0f },    def::math::Vector3{ 0.0f, 1.0f, 0.0f, 1.0f },    def::math::Vector3{ 0.0f, 0.0f, 0.0f, 1.0f },		def::math::Tex3{ 0.0f, 1.0f, 1.0f },		def::math::Tex3{ 1.0f, 0.0f, 1.0f },		def::math::Tex3{ 1.0f, 1.0f, 1.0f } },

			// TOP             																			   
			def::gfx3d::Triangle3f{ def::math::Vector3{ 0.0f, 1.0f, 0.0f, 1.0f },    def::math::Vector3{ 0.0f, 1.0f, 1.0f, 1.0f },    def::math::Vector3{ 1.0f, 1.0f, 1.0f, 1.0f },		def::math::Tex3{ 0.0f, 1.0f, 1.0f },		def::math::Tex3{ 0.0f, 0.0f, 1.0f },		def::math::Tex3{ 1.0f, 0.0f, 1.0f } },
			def::gfx3d::Triangle3f{ def::math::Vector3{ 0.0f, 1.0f, 0.0f, 1.0f },    def::math::Vector3{ 1.0f, 1.0f, 1.0f, 1.0f },    def::math::Vector3{ 1.0f, 1.0f, 0.0f, 1.0f },		def::math::Tex3{ 0.0f, 1.0f, 1.0f },		def::math::Tex3{ 1.0f, 0.0f, 1.0f },		def::math::Tex3{ 1.0f, 1.0f, 1.0f } },

			// BOTTOM          																			  
			def::gfx3d::Triangle3f{ def::math::Vector3{ 1.0f, 0.0f, 1.0f, 1.0f },    def::math::Vector3{ 0.0f, 0.0f, 1.0f, 1.0f },    def::math::Vector3{ 0.0f, 0.0f, 0.0f, 1.0f },		def::math::Tex3{ 0.0f, 1.0f, 1.0f },		def::math::Tex3{ 0.0f, 0.0f, 1.0f },		def::math::Tex3{ 1.0f, 0.0f, 1.0f } },
			def::gfx3d::Triangle3f{ def::math::Vector3{ 1.0f, 0.0f, 1.0f, 1.0f },    def::math::Vector3{ 0.0f, 0.0f, 0.0f, 1.0f },    def::math::Vector3{ 1.0f, 0.0f, 0.0f, 1.0f },		def::math::Tex3{ 0.0f, 1.0f, 1.0f },		def::math::Tex3{ 1.0f, 0.0f, 1.0f },		def::math::Tex3{ 1.0f, 1.0f, 1.0f } },
#endif
		};

		return m;
	}

	Engine::Engine(
		const def::vi2d& viewPort,
		float fovRad,
		float aspectRatio,
		float near,
		float far
	) : m_ViewPort(viewPort)
	{
		m_DepthBuffer = new float[viewPort.x * viewPort.y];
		m_Proj.projection(fovRad, aspectRatio, near, far);

		m_DrawTarget = new def::Graphic(viewPort.x, viewPort.y);

		m_WorldMat.identity();
	}

	Engine::~Engine()
	{
		delete[] m_DepthBuffer;
		delete m_TargetTexture;
		delete m_DrawTarget;
	}

	def::Graphic* Engine::DrawTarget() { return m_DrawTarget; }
	math::Matrix& Engine::WorldMatrix() { return m_WorldMat; }

	Sprite::SampleMethod& Engine::SampleMethod() { return m_TexSampleMethod; }
	Sprite::WrapMethod& Engine::WrapMethod() { return m_TexWrapMethod; }

	void Engine::Begin()
	{
		m_Model.tris.push_back({});
		m_FormingTriangle = &m_Model.tris.back();
		m_Vertex = 0;
		m_TexCoord = 0;

		m_FormingTriangle->col[0] = m_TargetCol[0];
		m_FormingTriangle->col[1] = m_TargetCol[1];
		m_FormingTriangle->col[2] = m_TargetCol[2];
	}

	void Engine::End()
	{
		if (m_FormingTriangle)
		{
			m_Vertex = 0;
			m_TexCoord = 0;
			m_FormingTriangle = nullptr;
		}
	}

	void Engine::Texture(const def::Sprite* texture) { m_TargetTexture = texture; }
	void Engine::Color(const def::Pixel& p1, const def::Pixel& p2, const def::Pixel& p3) 
	{
		m_TargetCol[0] = p1;
		m_TargetCol[1] = p2;
		m_TargetCol[2] = p3;
	}

	void Engine::Vertex(const math::Vector3f& v)
	{
		if (m_FormingTriangle && m_Vertex < 3)
			m_FormingTriangle->p[m_Vertex++] = v;
	}

	void Engine::TexCoord(const math::Tex3f& t)
	{
		if (m_FormingTriangle && m_TexCoord < 3)
			m_FormingTriangle->t[m_TexCoord++] = t;
	}

	void Engine::SetModel(const Mesh& model)
	{
		m_Model.tris = model.tris;
	}

	void Engine::Clear(const def::Pixel& col, bool clearDrawTarget, bool clearDepthBuffer)
	{
		if (clearDrawTarget)
		{
			m_DrawTarget->sprite->SetPixelData(col);
			m_DrawTarget->UpdateTexture();
		}

		if (clearDepthBuffer)
		{
			for (int i = 0; i < m_ViewPort.x * m_ViewPort.y; i++)
				m_DepthBuffer[i] = 0.0f;
		}
	}

	void Engine::Update(Camera3f& camera, const math::Vector3f& up, const math::Vector3f& target)
	{
		math::Matrix viewMat;
		viewMat.view(up, target, camera.pos, camera.dir, camera.rot);

		m_TrianglesToRaster.clear();

		for (const auto& tri : m_Model.tris)
		{
			Triangle3f transformed, viewed;

			for (size_t i = 0; i < 3; i++)
			{
				transformed.p[i] = m_WorldMat * tri.p[i];
				transformed.t[i] = tri.t[i];
				transformed.col[i] = m_TargetCol[i];
			}

			math::Vector3f line1 = transformed.p[1] - transformed.p[0];
			math::Vector3f line2 = transformed.p[2] - transformed.p[0];

			math::Vector3f normal = line1.cross(line2).norm();

			math::Vector3f ray = transformed.p[0] - camera.pos;

			if (normal.dot(ray) < 0.0f)
			{
				math::Vector3f lightDir = math::Vector3f(0.0f, 1.0f, -1.0f).norm();

				float dp = std::max(0.1f, lightDir.dot(normal));

				for (size_t i = 0; i < 3; i++)
				{
					viewed.p[i] = viewMat * transformed.p[i];
					viewed.t[i] = transformed.t[i];
					viewed.col[i].r = transformed.col[i].r * dp;
					viewed.col[i].g = transformed.col[i].g * dp;
					viewed.col[i].b = transformed.col[i].b * dp;
				}

				Triangle3f clipped[2];
				uint32_t clippedCount = viewed.clip_against_plane({ 0.0f, 0.0f, 0.1f }, { 0.0f, 0.0f, 1.0f }, clipped[0], clipped[1]);

				for (uint32_t n = 0; n < clippedCount; n++)
				{
					Triangle3f projected;

					for (size_t i = 0; i < 3; i++)
					{
						projected.p[i] = m_Proj * clipped[n].p[i];
						projected.t[i] = clipped[n].t[i];
						projected.col[i] = clipped[n].col[i];

						projected.p[i] /= projected.p[i].w;

						projected.t[i].u /= projected.p[i].w;
						projected.t[i].v /= projected.p[i].w;
						projected.t[i].w = 1.0f / projected.p[i].w;

						projected.p[i].x = (-projected.p[i].x + 1.0f) * 0.5f * (float)m_ViewPort.x;
						projected.p[i].y = (-projected.p[i].y + 1.0f) * 0.5f * (float)m_ViewPort.y;
					}

					m_TrianglesToRaster.push_back(projected);
				}
			}
		}

		m_TargetCol[0] = def::WHITE;
		m_TargetCol[1] = def::WHITE;
		m_TargetCol[2] = def::WHITE;

		m_Model.tris.clear();
	}

	void Engine::Draw()
	{
		for (auto& tri : m_TrianglesToRaster)
		{
			Triangle3f clipped[2];

			std::list<Triangle3f> tris = { tri };
			uint32_t newTris = 1;

			for (int p = 0; p < 4; p++)
			{
				uint32_t trisToAdd = 0;

				while (newTris > 0)
				{
					Triangle3f t = tris.front();
					tris.pop_front();
					newTris--;

					switch (p)
					{
					case 0:	trisToAdd = t.clip_against_plane({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, clipped[0], clipped[1]); break;
					case 1:	trisToAdd = t.clip_against_plane({ 0.0f, (float)m_ViewPort.y - 1.0f, 0.0f }, { 0.0f, -1.0f, 0.0f }, clipped[0], clipped[1]); break;
					case 2:	trisToAdd = t.clip_against_plane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, clipped[0], clipped[1]); break;
					case 3:	trisToAdd = t.clip_against_plane({ (float)m_ViewPort.x - 1.0f, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f }, clipped[0], clipped[1]); break;
					}

					for (int n = 0; n < trisToAdd; n++)
						tris.push_back(clipped[n]);
				}

				newTris = tris.size();
			}

			for (const auto& t : tris)
			{
				DrawTexturedTriangle(
					{ (int)t.p[0].x, (int)t.p[0].y }, t.t[0],
					{ (int)t.p[1].x, (int)t.p[1].y }, t.t[1],
					{ (int)t.p[2].x, (int)t.p[2].y }, t.t[2],
					t.col[0],
					t.col[1],
					t.col[2],
					m_TargetTexture
				);
			}
		}

		m_TargetTexture = nullptr;
		m_DrawTarget->UpdateTexture();
	}

	void Engine::DrawTexturedTriangle(
		def::vi2d p1, math::Tex3f t1,
		def::vi2d p2, math::Tex3f t2,
		def::vi2d p3, math::Tex3f t3,
		const def::Pixel& col1,
		const def::Pixel& col2,
		const def::Pixel& col3,
		const def::Sprite* spr
	)
	{
		if (p2.y < p1.y)
		{
			p1.swap(p2);
			t1.swap(t2);
		}

		if (p3.y < p1.y)
		{
			p1.swap(p3);
			t1.swap(t3);
		}

		if (p3.y < p2.y)
		{
			p2.swap(p3);
			t2.swap(t3);
		}

		def::vi2d dp1 = p2 - p1;
		def::vi2d dp2 = p3 - p1;

		math::Tex3f dt1 = t2 - t1;
		math::Tex3f dt2 = t3 - t1;

		math::Tex3f tex = { 0.0f, 0.0f, 0.0f };

		float dax_step = 0.0f, dbx_step = 0.0f;
		math::Tex3f dt1_step, dt2_step;
		dt1_step.w = 0.0f; dt2_step.w = 0.0f;

		if (dp2.y > 0.0f) dbx_step = dp2.x / (float)abs(dp2.y);
		if (dp2.y > 0.0f) dt2_step = dt2 / (float)abs(dp2.y);

		if (dp1.y > 0.0f)
		{
			dax_step = dp1.x / (float)abs(dp1.y);
			dt1_step = dt1 / (float)abs(dp1.y);

			for (int i = p1.y; i <= p2.y; i++)
			{
				int ax = p1.x + float(i - p1.y) * dax_step;
				int bx = p1.x + float(i - p1.y) * dbx_step;

				math::Tex3f texStart = t1 + float(i - p1.y) * dt1_step;
				math::Tex3f texEnd = t1 + float(i - p1.y) * dt2_step;

				if (ax > bx)
				{
					std::swap(ax, bx);
					texStart.swap(texEnd);
				}

				tex = texStart;

				float tstep = 1.0f / float(bx - ax);
				float t = 0.0f;

				for (int j = ax; j < bx; j++)
				{
					tex = (1.0f - t) * texStart + t * texEnd;

					if (tex.w > m_DepthBuffer[i * m_ViewPort.x + j])
					{
						m_DrawTarget->sprite->SetPixel(
							{ j, i },
							spr ?
								spr->Sample({ tex.u / tex.w, tex.v / tex.w }, m_TexSampleMethod, m_TexWrapMethod) :
								col1.mix(col2, tex.u / tex.w).mix(col3, tex.v / tex.w)
						);

						m_DepthBuffer[i * m_ViewPort.x + j] = tex.w;
					}

					t += tstep;
				}

			}
		}

		dp1 = p3 - p2;
		dt1 = t3 - t2;

		if (dp2.y > 0) dbx_step = dp2.x / (float)abs(dp2.y);

		dt1_step.u = 0.0f;
		dt1_step.v = 0.0f;

		if (dp1.y > 0.0f)
		{
			dax_step = dp1.x / (float)abs(dp1.y);
			dt1_step = dt1 / (float)abs(dp1.y);

			for (int i = p2.y; i <= p3.y; i++)
			{
				int ax = p2.x + float(i - p2.y) * dax_step;
				int bx = p1.x + float(i - p1.y) * dbx_step;

				math::Tex3f texStart = t2 + float(i - p2.y) * dt1_step;
				math::Tex3f texEnd = t1 + float(i - p1.y) * dt2_step;

				if (ax > bx)
				{
					std::swap(ax, bx);
					texStart.swap(texEnd);
				}

				tex = texStart;

				float tstep = 1.0f / float(bx - ax);
				float t = 0.0f;

				for (int j = ax; j < bx; j++)
				{
					tex = (1.0f - t) * texStart + t * texEnd;

					if (tex.w > m_DepthBuffer[i * m_ViewPort.x + j])
					{
						m_DrawTarget->sprite->SetPixel(
							{ j, i },
							spr ?
								spr->Sample({ tex.u / tex.w, tex.v / tex.w }, m_TexSampleMethod, m_TexWrapMethod) :
								col1.mix(col2, tex.u / tex.w).mix(col3, tex.v / tex.w)
						);

						m_DepthBuffer[i * m_ViewPort.x + j] = tex.w;
					}

					t += tstep;
				}
			}
		}
	}

#endif
}
