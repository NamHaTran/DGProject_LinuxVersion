#include "CommandCheck.h"
#include <string>

namespace postProcessKey
{
	std::string exitKey1("Exit"), exitKey2("exit"), exitKey3("EXIT"), exitKey4("End"), exitKey5("end"), exitKey6("END");
	
	/*Function return true if Exit is available*/
	bool checkExit(std::string cmd)
	{
		bool trigger(false);
		if ((cmd.compare(exitKey1) == 0) || (cmd.compare(exitKey2) == 0) || (cmd.compare(exitKey3) == 0) || (cmd.compare(exitKey4) == 0) || (cmd.compare(exitKey5) == 0) || (cmd.compare(exitKey6) == 0))
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
}

namespace processKey
{
	std::string runningKey1("DG2D"), runningKey2("dg2d");

	/*Function return true if DG2D is available*/
	bool checkDGRun(std::string cmd)
	{
		bool trigger(false);
		if ((cmd.compare(runningKey1) == 0) || (cmd.compare(runningKey2) == 0))
		{
			trigger = true;
		}
		return trigger;
	}
}
