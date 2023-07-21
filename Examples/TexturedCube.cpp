/*
*	BSD 3-Clause License

	Copyright (c) 2023, Alex

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are met:

	1. Redistributions of source code must retain the above copyright notice, this
	   list of conditions and the following disclaimer.

	2. Redistributions in binary form must reproduce the above copyright notice,
	   this list of conditions and the following disclaimer in the documentation
	   and/or other materials provided with the distribution.

	3. Neither the name of the copyright holder nor the names of its
	   contributors may be used to endorse or promote products derived from
	   this software without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
	FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
	DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
	SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
	CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
	OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
	OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#define DGE_APPLICATION
#include "../defGameEngine.h"

#define DEF_ENGINE_3D
#include "../defEngine3D.hpp"

class Engine3D : public def::GameEngine
{
public:
	Engine3D()
	{
		SetTitle("3D Engine");
	}

private:
	float speed = 8.0f;
	std::unique_ptr<def::gfx3d::Engine> engine3d;

	def::gfx3d::Camera3f camera;

	def::gfx3d::Mesh cube;
	std::unique_ptr<def::Sprite> tex;

	def::vu2d viewport;

protected:
	bool OnUserCreate() override
	{
		viewport = ScreenSize();

		engine3d = std::make_unique<def::gfx3d::Engine>(viewport, 3.14159f / 2.0f, (float)viewport.x / (float)viewport.y, 0.001f, 1000.0f);

		tex = std::make_unique<def::Sprite>("Assets/eye.png");

		cube = def::gfx3d::Mesh::MakeCube();

		return true;
	}

	float theta = 0.0f;

	bool OnUserUpdate(float dt) override
	{
		if (GetKey(def::Key::UP).held) camera.move_by(def::math::Vector3f(0.0f, speed * dt));
		if (GetKey(def::Key::DOWN).held) camera.move_by(def::math::Vector3f(0.0f, -speed * dt));

		def::math::Vector3f forward = camera.lookDir * speed * dt;

		if (GetKey(def::Key::W).held) camera.move_by(forward);
		if (GetKey(def::Key::S).held) camera.move_by(forward * -1.0f);
		if (GetKey(def::Key::A).held) camera.rotate_by(speed * 0.5f * dt);
		if (GetKey(def::Key::D).held) camera.rotate_by(-speed * 0.5f * dt);

		engine3d->Texture(tex.get());
		engine3d->SetModel(cube);

		def::math::Matrix rotMatZ, rotMatX, rotMatXZ;

		rotMatZ.rotate_z(theta);
		rotMatX.rotate_x(theta);

		rotMatXZ = rotMatZ * rotMatX;
		engine3d->WorldMatrix() = rotMatXZ;

		engine3d->Update(camera);

		engine3d->Clear(def::CYAN);
		engine3d->Draw();

		SetDrawTarget(engine3d->DrawTarget());

		theta += dt;

		return true;
	}
};

int main()
{
	Engine3D demo;

	demo.Construct(256, 240, 4, 4);
	demo.Run();

	return 0;
}
