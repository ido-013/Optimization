#include "Profiler.h"

#include <iostream>

MyProfiler::Block::Block(const std::string& _name, Block* _parent) : name(_name), parent(_parent)
{
	start = std::chrono::steady_clock::now();
}

MyProfiler::Block::~Block()
{
	//delete all the children
	for (auto it : children)
		delete it;

	children.clear();
}

void MyProfiler::Block::End()
{
	end = std::chrono::steady_clock::now();
}

double MyProfiler::Block::GetSeconds() const
{
	return std::chrono::duration<double>(end - start).count();
}

MyProfiler::Block* MyProfiler::Block::AddChild(const std::string& _name)
{
	Block* nBlock = new Block(_name, this);
	children.push_back(nBlock);
	return nBlock;
}

void MyProfiler::Block::Dump(int _n) const
{
	//print correct ammount of tabs
	for (int i = 0; i < _n; i++)
	{
		std::cout << "\t";
	}

	//print name and seconds
	std::cout << name << " in " << GetSeconds() << " seconds" << std::endl;

	//print children
	for (const auto* c : children)
	{
		c->Dump(_n + 1);
	}
}

MyProfiler::Profiler* MyProfiler::Profiler::ptr = nullptr;

MyProfiler::Profiler::~Profiler()
{
	Clear();
}

void MyProfiler::Profiler::StartBlock(std::string _name)
{
	if (!current)
		current = new Block(_name); //I am root
	else
		current = current->AddChild(_name);
}

void MyProfiler::Profiler::End()
{
	//stop counting time on current block
	current->End();

	//go 1 step back
	Block* parent = current->GetParent();

	//if no parent. Push current to fullyFinished
	if (!parent)
		fullyFinishedBlocks.push_back(current);

	current = parent;
}

void MyProfiler::Profiler::Dump()
{
	for (const auto* b : fullyFinishedBlocks)
	{
		b->Dump();
	}
}

void MyProfiler::Profiler::Clear()
{
	//iterate end() UNTIL current in nullptr
	while (current)
		End();

	//delete all the finished nodes
	for (auto it : fullyFinishedBlocks)
		delete it;

	fullyFinishedBlocks.clear();
}

MyProfiler::Profiler* MyProfiler::Profiler::GetPtr()
{
	if (!ptr)
		ptr = new Profiler;

	return ptr;
}

void MyProfiler::Profiler::DeletePtr()
{
	if (ptr)
		delete ptr;

	ptr = nullptr;
}
