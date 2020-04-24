#include "CommandCheck.h"
#include <string>

namespace postProcessKey
{
	/*Function return true if Exit is available*/
	bool checkExit(std::string cmd)
	{
		bool trigger(false);
        if ((cmd.compare("Exit") == 0) || (cmd.compare("exit") == 0) || (cmd.compare("EXIT") == 0) || (cmd.compare("End") == 0) || (cmd.compare("end") == 0) || (cmd.compare("END") == 0))
		{
			trigger = true;
		}
		return trigger;
	}
}

namespace preProcessKey
{
	/*Function return true if UnvToDG is available*/
	bool checkUnvReader(std::string cmd)
	{
		bool trigger(false);
		if ((cmd.compare("UnvToDG") == 0) || (cmd.compare("unvtodg") == 0) || (cmd.compare("UNVTODG") == 0))
		{
			trigger = true;
		}
		return trigger;
	}

	/*Function return true if UnvReaderHelp is available*/
	bool checkUnvHelper(std::string cmd)
	{
		bool trigger(false);
		if ((cmd.compare("UnvReaderHelp") == 0) || (cmd.compare("unvreaderhelp") == 0) || (cmd.compare("unvreaderhelp") == 0))
		{
			trigger = true;
		}
		return trigger;
	}
	/*Function return true if bcHelp is available*/
	bool checkBCsHelper(std::string cmd)
	{
		bool trigger(false);
		if ((cmd.compare("bchelp") == 0) || (cmd.compare("BcHelp") == 0))
		{
			trigger = true;
		}
		return trigger;
	}

	bool reSubmit(std::string cmd)
	{
		bool trigger(false);
		if ((cmd.compare("resubmit") == 0) || (cmd.compare("reSubmit") == 0))
		{
			trigger = true;
		}
		return trigger;
	}

	bool mappResults(std::string cmd)
	{
		bool trigger(false);
        if ((cmd.compare("mapresults") == 0) || (cmd.compare("mapResults") == 0))
		{
			trigger = true;
		}
		return trigger;
	}

    bool exportMeshToMetis(std::string cmd)
    {
        bool trigger(false);
        if ((cmd.compare("DGMeshToMetis") == 0) || (cmd.compare("dgmeshtometis") == 0))
        {
            trigger = true;
        }
        return trigger;
    }

    bool testMeshPartitionResult(std::string cmd)
    {
        bool trigger(false);
        if ((cmd.compare("testmeshpartitionresult") == 0) || (cmd.compare("testMeshPartitionResult") == 0))
        {
            trigger = true;
        }
        return trigger;
    }

    bool decomposeCase(std::string cmd)
    {
        bool trigger(false);
        if ((cmd.compare("decomposecase") == 0) || (cmd.compare("decomposeCase") == 0))
        {
            trigger = true;
        }
        return trigger;
    }

    bool reconstructLatestTime(std::string cmd)
    {
        bool trigger(false);
        if ((cmd.compare("reconstructCase") == 0) || (cmd.compare("reconstructcase") == 0))
        {
            trigger = true;
        }
        return trigger;
    }

    bool checkPartitionedMesh(std::string cmd)
    {
        bool trigger(false);
        if ((cmd.compare("checkParMesh") == 0) || (cmd.compare("checkparmesh") == 0))
        {
            trigger = true;
        }
        return trigger;
    }

	namespace debug
	{
		bool checkElement(std::string cmd)
		{
			bool trigger(false);
			if ((cmd.compare("checkElementInfor") == 0) || (cmd.compare("checkelementinfor") == 0))
			{
				trigger = true;
			}
			return trigger;
		}
	}
}

namespace processKey
{
	/*Function return true if DG2D is available*/
	bool checkDGRun(std::string cmd)
	{
		bool trigger(false);
        if ((cmd.compare("DG2D") == 0) || (cmd.compare("dg2d") == 0))
		{
			trigger = true;
		}
		return trigger;
	}
}
