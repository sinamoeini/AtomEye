/******************************************************************************
 * The code to dynamically load HAPI.DLL													*
 * You may compile and link it to your program											*
 * or source it using '#include "hcload.c"'												*
 ******************************************************************************/

#include "windows.h"
#define HC_LOAD_CODE 
#include "hc.h"
#undef HC_LOAD_CODE

#ifndef TRUE
	#define TRUE 1
	#define FALSE 0
#endif


BOOL _hcLoadError(char* fnName, char* libName)
{
#define _hcMaxMess 180
char sz[_hcMaxMess];

if (lstrlen(libName) > 80) 
     wsprintf(sz,"Error during loading hcAPI: Missing function %s in library",fnName);
else
     wsprintf(sz,"Error during loading hcAPI: Missing function %s in library %s",fnName,libName);

MessageBox(NULL,sz,"hcAPI Load",MB_OK);

return FALSE;

#undef _hcMaxMess
}


BOOL LoadHAPI(LPSTR szN)
	 {
	  HINSTANCE hinst;

	  BOOL res=TRUE;

	  hinst=LoadLibrary(szN);
	  if (hinst==NULL) 
	  { MessageBox(NULL,"error while loading HAPI\n","HAPI Init",MB_OK); 
	    return FALSE;
	   };



	 hcInitAPI=(T_hcInitAPI*)GetProcAddress(hinst,"hcInitAPI");
      if (hcInitAPI == NULL) res=_hcLoadError("hcInitAPI",szN);
		else (hcInitAPI)();


	 hcConnect=(T_hcConnect*)GetProcAddress(hinst,"hcConnect");
      if ( hcConnect == NULL) res=_hcLoadError("hcConnect",szN);


	 hcDisconnect=(T_hcDisconnect*)GetProcAddress(hinst,"hcDisconnect");
      if ( hcDisconnect == NULL) res=_hcLoadError("hcDisconnect",szN);

	 hcExit=(T_hcExit*)GetProcAddress(hinst,"hcExit");
      if ( hcExit == NULL) res=_hcLoadError("hcExit",szN);

	 hcExecTxt=(T_hcExecTxt*)GetProcAddress(hinst,"hcExecTxt");
      if ( hcExecTxt == NULL) res=_hcLoadError("hcExecTxt",szN);

	 hcExecBin=(T_hcExecBin*)GetProcAddress(hinst,"hcExecBin");
      if ( hcExecBin == NULL) res=_hcLoadError("hcExecBin",szN);

	 hcQueryTxt=(T_hcQueryTxt*)GetProcAddress(hinst,"hcQueryTxt");
      if ( hcQueryTxt == NULL) res=_hcLoadError("hcQueryTxt",szN);

	 hcQueryBin=(T_hcQueryBin*)GetProcAddress(hinst,"hcQueryBin");
      if ( hcQueryBin == NULL) res=_hcLoadError("hcQueryBin",szN);

 	 hcGetInt=(T_hcGetInt*)GetProcAddress(hinst,"hcGetInt");
      if ( hcGetInt == NULL) res=_hcLoadError("hcGetInt",szN);

	 hcGetReal=(T_hcGetReal*)GetProcAddress(hinst,"hcGetReal");
      if ( hcGetReal == NULL) res=_hcLoadError("hcGetReal",szN);

	 hcGetIntVec=(T_hcGetIntVec*)GetProcAddress(hinst,"hcGetIntVec");
      if (  hcGetIntVec == NULL) res=_hcLoadError("hcGetIntVec",szN);

	 hcGetIntArr=(T_hcGetIntArr*)GetProcAddress(hinst,"hcGetIntArr");
      if (  hcGetIntArr == NULL) res=_hcLoadError("hcGetIntArr",szN);


	 hcGetRealVec=(T_hcGetRealVec*)GetProcAddress(hinst,"hcGetRealVec");
      if (  hcGetRealVec == NULL) res=_hcLoadError("hcGetRealVec",szN);


	 hcGetRealArr=(T_hcGetRealArr*)GetProcAddress(hinst,"hcGetRealArr");
      if (  hcGetRealArr == NULL) res=_hcLoadError("hcGetRealArr",szN);


	 hcGetIntVecElm=(T_hcGetIntVecElm*)GetProcAddress(hinst,"hcGetIntVecElm");
      if (  hcGetIntVecElm == NULL) res=_hcLoadError("hcGetIntVecElm",szN);


	 hcGetRealVecElm=(T_hcGetRealVecElm*)GetProcAddress(hinst,"hcGetRealVecElm");
      if (  hcGetRealVecElm == NULL) res=_hcLoadError("hcGetRealVecElm",szN);


	 hcGetIntArrElm=(T_hcGetIntArrElm*)GetProcAddress(hinst,"hcGetIntArrElm");
      if (  hcGetIntArrElm == NULL) res=_hcLoadError("hcGetIntArrElm",szN);


	 hcGetRealArrElm=(T_hcGetRealArrElm*)GetProcAddress(hinst,"hcGetRealArrElm");
      if (  hcGetRealArrElm == NULL) res=_hcLoadError("hcGetRealArrElm",szN);


	 hcGetRealVecXYZ=(T_hcGetRealVecXYZ*)GetProcAddress(hinst,"hcGetRealVecXYZ");
      if (  hcGetRealVecXYZ == NULL) res=_hcLoadError("hcGetRealVecXYZ",szN);


	 hcGetRealArrXYZ=(T_hcGetRealArrXYZ*)GetProcAddress(hinst,"hcGetRealArrXYZ");
      if (  hcGetRealArrXYZ == NULL) res=_hcLoadError("hcGetRealArrXYZ",szN);


	 hcGetStr=(T_hcGetStr*)GetProcAddress(hinst,"hcGetStr");
      if (  hcGetStr == NULL) res=_hcLoadError("hcGetStr",szN);


	 hcGetStrVecElm=(T_hcGetStrVecElm*)GetProcAddress(hinst,"hcGetStrVecElm");
      if (  hcGetStrVecElm == NULL) res=_hcLoadError("hcGetStrVecElm",szN);


	 hcGetStrArrElm=(T_hcGetStrArrElm*)GetProcAddress(hinst,"hcGetStrArrElm");
      if (  hcGetStrArrElm == NULL) res=_hcLoadError("hcGetStrArrElm",szN);


	 hcGetBlock=(T_hcGetBlock*)GetProcAddress(hinst,"hcGetBlock");
      if (  hcGetBlock == NULL) res=_hcLoadError("hcGetBlock",szN);


	 hcSetInt=(T_hcSetInt*)GetProcAddress(hinst,"hcSetInt");
      if (  hcSetInt == NULL) res=_hcLoadError("hcSetInt",szN);


	 hcSetReal=(T_hcSetReal*)GetProcAddress(hinst,"hcSetReal");
      if (  hcSetReal == NULL) res=_hcLoadError("hcSetReal",szN);


	 hcSetIntVec=(T_hcSetIntVec*)GetProcAddress(hinst,"hcSetIntVec");
      if (  hcSetIntVec == NULL) res=_hcLoadError("hcSetIntVec",szN);


	 hcSetIntArr=(T_hcSetIntArr*)GetProcAddress(hinst,"hcSetIntArr");
      if (  hcSetIntArr == NULL) res=_hcLoadError("hcSetIntArr",szN);


	 hcSetRealVec=(T_hcSetRealVec*)GetProcAddress(hinst,"hcSetRealVec");
      if (  hcSetRealVec == NULL) res=_hcLoadError("hcSetRealVec",szN);


	 hcSetRealArr=(T_hcSetRealArr*)GetProcAddress(hinst,"hcSetRealArr");
      if (  hcSetRealArr == NULL) res=_hcLoadError("hcSetRealArr",szN);


	 hcSetIntVecElm=(T_hcSetIntVecElm*)GetProcAddress(hinst,"hcSetIntVecElm");
      if (  hcSetIntVecElm == NULL) res=_hcLoadError("hcSetIntVecElm",szN);


	 hcSetRealVecElm=(T_hcSetRealVecElm*)GetProcAddress(hinst,"hcSetRealVecElm");
      if (  hcSetRealVecElm == NULL) res=_hcLoadError("hcSetRealVecElm",szN);


	 hcSetIntArrElm=(T_hcSetIntArrElm*)GetProcAddress(hinst,"hcSetIntArrElm");
      if (  hcSetIntArrElm == NULL) res=_hcLoadError("hcSetIntArrElm",szN);


	 hcSetRealArrElm=(T_hcSetRealArrElm*)GetProcAddress(hinst,"hcSetRealArrElm");
      if (  hcSetRealArrElm == NULL) res=_hcLoadError("hcSetRealArrElm",szN);


	 hcSetRealArrXYZ=(T_hcSetRealArrXYZ*)GetProcAddress(hinst,"hcSetRealArrXYZ");
      if (  hcSetRealArrXYZ == NULL) res=_hcLoadError("hcSetRealArrXYZ",szN);


	 hcSetRealVecXYZ=(T_hcSetRealVecXYZ*)GetProcAddress(hinst,"hcSetRealVecXYZ");
      if (  hcSetRealVecXYZ == NULL) res=_hcLoadError("hcSetRealVecXYZ",szN);


	 hcSetStr=(T_hcSetStr*)GetProcAddress(hinst,"hcSetStr");
      if (  hcSetStr == NULL) res=_hcLoadError("hcSetStr",szN);


	 hcSetStrVecElm=(T_hcSetStrVecElm*)GetProcAddress(hinst,"hcSetStrVecElm");
      if (  hcSetStrVecElm == NULL) res=_hcLoadError("hcSetStrVecElm",szN);


	 hcSetStrArrElm=(T_hcSetStrArrElm*)GetProcAddress(hinst,"hcSetStrArrElm");
      if (  hcSetStrArrElm == NULL) res=_hcLoadError("hcSetStrArrElm",szN);

	 hcSetBlock=(T_hcSetBlock*)GetProcAddress(hinst,"hcSetBlock");
      if (  hcSetBlock == NULL) res=_hcLoadError("hcSetBlock",szN);

     hcAlloc=(T_hcAlloc*)GetProcAddress(hinst,"hcAlloc");
      if ( hcAlloc== NULL) res=_hcLoadError("hcAlloc",szN);

	 hcFree=(T_hcFree*)GetProcAddress(hinst,"hcFree");
      if ( hcFree == NULL) res=_hcLoadError("hcFree",szN);

	 hcShowMessage=(T_hcShowMessage*)GetProcAddress(hinst,"hcShowMessage");
      if ( hcShowMessage == NULL) res=_hcLoadError("hcShowMessage",szN);

	 hcLastError=(T_hcLastError*)GetProcAddress(hinst,"hcLastError");
      if ( hcLastError == NULL) res=_hcLoadError("hcLastError",szN);

	 hcSetTimeouts=(T_hcSetTimeouts*)GetProcAddress(hinst,"hcSetTimeouts");
      if ( hcSetTimeouts == NULL) res=_hcLoadError("hcSetTimeouts",szN);

 	 hcGetErrorAction=(T_hcGetErrorAction*)GetProcAddress(hinst,"hcGetErrorAction");
      if ( hcGetErrorAction == NULL) res=_hcLoadError("hcGetErrorAction",szN);
	 
	 hcSetErrorAction=(T_hcSetErrorAction*)GetProcAddress(hinst,"hcSetErrorAction");
      if ( hcSetErrorAction == NULL) res=_hcLoadError("hcSetErrorAction",szN);

	 hcNotifyDataAvail=(T_hcNotifyDataAvail*)GetProcAddress(hinst,"hcNotifyDataAvail");
      if ( hcNotifyDataAvail == NULL) res=_hcLoadError("hcNotifyDataAvail",szN);

	 hcNotifyDataAvail=(T_hcNotifyDataAvail*)GetProcAddress(hinst,"hcNotifyDataAvail");
      if ( hcNotifyDataAvail == NULL) res=_hcLoadError("hcNotifyDataAvail",szN);

	 hcNotifyStart=(T_hcNotifyStart*)GetProcAddress(hinst,"hcNotifyStart");
      if ( hcNotifyStart == NULL) res=_hcLoadError("hcNotifyStart",szN);

	 hcNotifyStop=(T_hcNotifyStop*)GetProcAddress(hinst,"hcNotifyStop");
      if ( hcNotifyStop == NULL) res=_hcLoadError("hcNotifyStop",szN);

	 hcNotifySetup=(T_hcNotifySetup*)GetProcAddress(hinst,"hcNotifySetup");
      if ( hcNotifySetup == NULL) res=_hcLoadError("hcNotifySetup",szN);

	 hcGetNotifyData=(T_hcGetNotifyData*)GetProcAddress(hinst,"hcGetNotifyData");
      if ( hcGetNotifyData == NULL) res=_hcLoadError("hcGetNotifyData",szN);

return res;
	 }
