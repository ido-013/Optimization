#pragma once

#include <string>
#include <chrono>
#include <list>

namespace MyProfiler
{
	class Block
	{
	private:
		std::string name;
		std::chrono::steady_clock::time_point start;
		std::chrono::steady_clock::time_point end;

		std::list<Block*> children;
		Block* parent;

	public:
		Block(const std::string& _name, Block* _parent = nullptr);
		~Block();

		void End();

		double GetSeconds() const;
		std::string GetName() const { return name; }
		Block* GetParent() const { return parent; }

		Block* AddChild(const std::string& _name);

		void Dump(int _n = 0) const;
	};

	class Profiler
	{
	private:
		static Profiler* ptr;
		Profiler() = default;
		Profiler(const Profiler&) = default;
		~Profiler();
		Profiler& operator=(const Profiler&) = default;

		Block* current = nullptr;

		std::list<Block*> fullyFinishedBlocks;

	public:
		void StartBlock(std::string _name);
		void End();
		void Dump();

		void Clear();

		//Singleton
		static Profiler* GetPtr();
		static void DeletePtr();
	};
}

#ifndef __FUNCTION_NAME__
#ifdef WIN32
#define __FUNCTION_NAME__ __FUNCTION__
#else
#define __FUNCTION_NAME__ __func__
#endif // Win32
#endif // !__FUNCTION_NAME__

#ifdef _DEBUG
#define DEBUG_PROFILER_START(x) MyProfiler::Profiler::GetPtr()->StartBlock(x)
#define DEBUG_PROFILER_END MyProfiler::Profiler::GetPtr()->End()
#define DEBUG_PROFILER_DUMP MyProfiler::Profiler::GetPtr()->Dump()
#define DEBUG_PROFILER_DELETE MyProfiler::Profiler::DeletePtr()
#else
#define DEBUG_PROFILER_START(x) //
#define DEBUG_PROFILER_END		//
#define DEBUG_PROFILER_DUMP		//
#define DEBUG_PROFILER_DELETE	//
#endif