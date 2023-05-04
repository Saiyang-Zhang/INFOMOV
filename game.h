// Template, IGAD version 3
// Get the latest version from: https://github.com/jbikker/tmpl8
// IGAD/NHTV/UU - Jacco Bikker - 2006-2023

#pragma once

namespace Tmpl8
{

class Game : public TheApp
{
public:
	void Init();
	__int64 Evaluate();
	void Shutdown();
	void Tick( float _DT );
	void MouseUp( int _Button ) { /* implement if you want to detect mouse button presses */ }
	void MouseDown( int _Button ) { /* implement if you want to detect mouse button presses */ }
	void MouseMove( int _X, int _Y ) { /* implement if you want to detect mouse movement */ }
	void KeyUp( int _Key ) { /* implement if you want to handle keys */ }
	void KeyDown( int _Key ) { /* implement if you want to handle keys */ }
};

}; // namespace Tmpl8