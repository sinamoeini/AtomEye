/****hc.h**********************************************************************
 * Header file for the CDK API												  *
 ******************************************************************************/

#ifndef HC_HEADER_INCLUDED

/* Basic type definitions */

typedef  BYTE*   LPB;
typedef  void*   LPV;
typedef  DWORD*  LPLW;
typedef  long    HSV;       //Type placeholder for the HSV - HyperChem State Variable
typedef	 unsigned long MID; //Type placeholder for the MID - Message IDentifier

/* Definitions */
#define		hcMaxNameSize _MAX_PATH
#define		hcMaxMessSize 256

/* Types of errors  (bits 0-3) */

#define errNO_ERROR			0		// no error

#define errFATAL			1       // FATAL error of unknown nature
#define errNON_FATAL		2	    // NONFATAL error of unknown nature

#define errIN_HCAPI			4       // FATAL error in HCAPI
#define errIN_SYS			0       // FATAL error while in system
									// currently not used


/* Actions defined when the error occurs */

#define errACTION_NO				 0 // do not perform any action on any error
#define errACTION_MESS_BOX			16 // display message box
#define errACTION_DISCONNECT		32 // disconnect
#define errACTION_EXIT 			    64 // exit
#define errDDE_REP					1  // report DDE errors 
#define errDDE_NO_REP				2  // do not report DDE errors

/* Type Definitions for user defined callback functions */

typedef VOID (*PFNB)(DWORD, char*, DWORD);  // function type for Binary Callback Functions
typedef VOID (*PFNX)(LPSTR, LPSTR);			// function type for Text Callback Functions

	
/* Functional types for functions */

typedef BOOL			(_stdcall T_hcInitAPI)(void);  
typedef BOOL 			(_stdcall T_hcConnect)(LPSTR lszCmd);
typedef BOOL			(_stdcall T_hcDisconnect)(void);
typedef	void			(_stdcall T_hcExit)(void);

typedef BOOL            (_stdcall T_hcExecTxt)(LPSTR script_cmd);
typedef BOOL			(_stdcall T_hcExecBin)(HSV cmd, LPV args, DWORD args_length);

typedef LPSTR			(_stdcall T_hcQueryTxt)(LPSTR var_name);
typedef LPV 			(_stdcall T_hcQueryBin)(HSV var,DWORD indx1,DWORD indx2,DWORD* resp_length);

typedef int				(_stdcall T_hcGetInt)(HSV var);
typedef double          (_stdcall T_hcGetReal)(HSV var);
typedef int				(_stdcall T_hcGetIntVec)(HSV var, int* buff, int max_length);
typedef int				(_stdcall T_hcGetIntArr)(HSV var, int* buff, int max_length);
typedef int             (_stdcall T_hcGetRealVec)(HSV var, double* buff,int max_length);
typedef int				(_stdcall T_hcGetRealArr)(HSV var, double* buff,int max_length);

typedef int				(_stdcall T_hcGetIntVecElm)(HSV var,int index);
typedef double          (_stdcall T_hcGetRealVecElm)(HSV var,int index);
typedef int				(_stdcall T_hcGetIntArrElm)(HSV var,int atom_index,int molecule_index);
typedef double          (_stdcall T_hcGetRealArrElm)(HSV var,int atom_index,int molecule_index);
typedef BOOL            (_stdcall T_hcGetRealArrXYZ)(HSV var,int atom_index,int molecule_index, double* x, double* y, double* z);
typedef BOOL            (_stdcall T_hcGetRealVecXYZ)(HSV var,int atom_index,double* x, double* y, double* z);
typedef int				(_stdcall T_hcGetStr)(HSV var,char* buff,int max_length);
typedef int				(_stdcall T_hcGetStrVecElm)(HSV var, int index, char* buff, int max_length);
typedef int				(_stdcall T_hcGetStrArrElm)(HSV var, int atom_index, int molecule_index, char* buff, int max_length);
typedef int				(_stdcall T_hcGetBlock)(HSV var, char* buff, int max_length);

typedef BOOL			(_stdcall T_hcSetInt)(HSV var,int value);
typedef	BOOL	        (_stdcall T_hcSetReal)(HSV var,double value);
typedef BOOL			(_stdcall T_hcSetIntVec)(HSV var, int* buff, int length);
typedef BOOL			(_stdcall T_hcSetIntArr)(HSV var, int* buff, int length);
typedef BOOL            (_stdcall T_hcSetRealVec)(HSV var, double* buff,int length);
typedef BOOL			(_stdcall T_hcSetRealArr)(HSV var, double* buff,int length);

typedef BOOL			(_stdcall T_hcSetIntVecElm)(HSV var,int index,int value);
typedef BOOL	        (_stdcall T_hcSetRealVecElm)(HSV var,int index,double value);
typedef BOOL			(_stdcall T_hcSetIntArrElm)(HSV var,int atom_index,int molecule_index,int value);
typedef BOOL		    (_stdcall T_hcSetRealArrElm)(HSV var,int atom_index,int molecule_index,double value);
typedef BOOL            (_stdcall T_hcSetRealArrXYZ)(HSV var,int atom_index,int molecule_index,double x, double y, double z);
typedef BOOL            (_stdcall T_hcSetRealVecXYZ)(HSV var,int index,double x, double y, double z);
typedef BOOL			(_stdcall T_hcSetStr)(HSV var,char* string);
typedef BOOL			(_stdcall T_hcSetStrVecElm)(HSV var, int index, char* string);
typedef BOOL			(_stdcall T_hcSetStrArrElm)(HSV var,int atom_index, int molecule_index, char* string);
typedef BOOL			(_stdcall T_hcSetBlock)(HSV var, char* buff, int length);

typedef BOOL            (_stdcall T_hcNotifyStart)(LPSTR var_name);
typedef BOOL            (_stdcall T_hcNotifyStop)(LPSTR var_name);
typedef BOOL			(_stdcall T_hcNotifySetup)(PFNB pCallBack,BOOL NotifyWithText);
typedef DWORD			(_stdcall T_hcNotifyDataAvail)(int *data_size);
typedef DWORD			(_stdcall T_hcGetNotifyData)(char* name,char *buffer, DWORD MaxBuffLength);

typedef void*           (_stdcall T_hcAlloc)(size_t n_bytes);
typedef void			(_stdcall T_hcFree)(void* pointer);
typedef void			(_stdcall T_hcShowMessage)(LPSTR message);
typedef void			(_stdcall T_hcSetTimeouts)(DWORD ExcTimeOut,DWORD QryTimeOut,DWORD RstTimeOut);

typedef int				(_stdcall T_hcLastError)(char* szLastErr);
typedef int				(_stdcall T_hcGetErrorAction)(void);
typedef void			(_stdcall T_hcSetErrorAction)(int err);

#ifdef HC_BUILD  /* hc.h is included during normal hcAPI build */

		T_hcInitAPI 	hcInitAPI;
		T_hcConnect 	hcConnect;
		T_hcDisconnect	hcDisconnect;
		T_hcExit 		hcExit;

		T_hcExecTxt 	hcExecTxt;
		T_hcExecBin 	hcExecBin;

		T_hcQueryTxt 	hcQueryTxt;
		T_hcQueryBin 	hcQueryBin;

		T_hcGetInt			hcGetInt;
		T_hcGetReal			hcGetReal;
		T_hcGetIntVec		hcGetIntVec;
		T_hcGetIntArr		hcGetIntArr;
		T_hcGetRealVec		hcGetRealVec;
		T_hcGetRealArr		hcGetRealArr;
		T_hcGetIntVecElm	hcGetIntVecElm;
		T_hcGetRealVecElm	hcGetRealVecElm;
		T_hcGetIntArrElm	hcGetIntArrElm;
		T_hcGetRealArrElm	hcGetRealArrElm;
		T_hcGetRealVecXYZ   hcGetRealVecXYZ;
		T_hcGetRealArrXYZ   hcGetRealArrXYZ;
		T_hcGetStr			hcGetStr;
		T_hcGetStrVecElm	hcGetStrVecElm;
		T_hcGetStrArrElm	hcGetStrArrElm;
		T_hcGetBlock		hcGetBlock;


		T_hcSetInt			hcSetInt;
		T_hcSetReal			hcSetReal;
		T_hcSetIntVec		hcSetIntVec;
		T_hcSetIntArr		hcSetIntArr;
		T_hcSetRealVec		hcSetRealVec;
		T_hcSetRealArr		hcSetRealArr;
		                   
		T_hcSetIntVecElm	hcSetIntVecElm;
		T_hcSetRealVecElm	hcSetRealVecElm;
		T_hcSetIntArrElm	hcSetIntArrElm;
		T_hcSetRealArrElm	hcSetRealArrElm;
		T_hcSetRealArrXYZ	hcSetRealArrXYZ;
		T_hcSetRealVecXYZ	hcSetRealVecXYZ;
		T_hcSetStr			hcSetStr;
		T_hcSetStrVecElm	hcSetStrVecElm;
		T_hcSetStrArrElm	hcSetStrArrElm;
		T_hcSetBlock		hcSetBlock;


		T_hcNotifyStart		hcNotifyStart;
		T_hcNotifyStop		hcNotifyStop;
		T_hcNotifySetup		hcNotifySetup;
		T_hcGetNotifyData	hcGetNotifyData;
		T_hcNotifyDataAvail hcNotifyDataAvail;
		
		T_hcAlloc 			hcAlloc;
		T_hcFree 			hcFree;
		T_hcShowMessage 	hcShowMessage;
		T_hcSetTimeouts     hcSetTimeouts;

		T_hcLastError 			hcLastError;
		T_hcGetErrorAction      hcGetErrorAction;
		T_hcSetErrorAction 		hcSetErrorAction;

#else   /* hc.h is included for user interface to hcAPI.DLL */

#ifdef HC_LOAD_CODE /* for hcload.c we do not need pointers */

#	define _hc_TYPE_PREFIX extern

#else 

#	define _hc_TYPE_PREFIX 

#endif
		BOOL LoadHAPI(LPSTR szN);

		_hc_TYPE_PREFIX T_hcInitAPI 	*hcInitAPI;
		_hc_TYPE_PREFIX T_hcConnect 	*hcConnect;
		_hc_TYPE_PREFIX T_hcDisconnect	*hcDisconnect;
		_hc_TYPE_PREFIX T_hcExit 		*hcExit;

		_hc_TYPE_PREFIX T_hcExecTxt 	*hcExecTxt;
		_hc_TYPE_PREFIX T_hcExecBin 	*hcExecBin;

		_hc_TYPE_PREFIX T_hcQueryTxt 	*hcQueryTxt;
		_hc_TYPE_PREFIX T_hcQueryBin 	*hcQueryBin;

		_hc_TYPE_PREFIX T_hcGetInt			*hcGetInt;
		_hc_TYPE_PREFIX T_hcGetReal			*hcGetReal;
		_hc_TYPE_PREFIX T_hcGetIntVec		*hcGetIntVec;
		_hc_TYPE_PREFIX T_hcGetIntArr		*hcGetIntArr;
		_hc_TYPE_PREFIX T_hcGetRealVec		*hcGetRealVec;
		_hc_TYPE_PREFIX T_hcGetRealArr		*hcGetRealArr;
		_hc_TYPE_PREFIX T_hcGetIntVecElm		*hcGetIntVecElm;
		_hc_TYPE_PREFIX T_hcGetRealVecElm	*hcGetRealVecElm;
		_hc_TYPE_PREFIX T_hcGetIntArrElm		*hcGetIntArrElm;
		_hc_TYPE_PREFIX T_hcGetRealArrElm	*hcGetRealArrElm;
		_hc_TYPE_PREFIX T_hcGetRealVecXYZ	*hcGetRealVecXYZ;
		_hc_TYPE_PREFIX T_hcGetRealArrXYZ	*hcGetRealArrXYZ;
		_hc_TYPE_PREFIX T_hcGetStr			*hcGetStr;
		_hc_TYPE_PREFIX T_hcGetStrVecElm		*hcGetStrVecElm;
		_hc_TYPE_PREFIX T_hcGetStrArrElm		*hcGetStrArrElm;
		_hc_TYPE_PREFIX T_hcGetBlock			*hcGetBlock;


		_hc_TYPE_PREFIX T_hcSetInt			*hcSetInt;
		_hc_TYPE_PREFIX T_hcSetReal			*hcSetReal;
		_hc_TYPE_PREFIX T_hcSetIntVec		*hcSetIntVec;
		_hc_TYPE_PREFIX T_hcSetIntArr		*hcSetIntArr;
		_hc_TYPE_PREFIX T_hcSetRealVec		*hcSetRealVec;
		_hc_TYPE_PREFIX T_hcSetRealArr		*hcSetRealArr;
		                 
		_hc_TYPE_PREFIX T_hcSetIntVecElm		*hcSetIntVecElm;
		_hc_TYPE_PREFIX T_hcSetRealVecElm	*hcSetRealVecElm;
		_hc_TYPE_PREFIX T_hcSetIntArrElm		*hcSetIntArrElm;
		_hc_TYPE_PREFIX T_hcSetRealArrElm	*hcSetRealArrElm;
		_hc_TYPE_PREFIX T_hcSetRealArrXYZ	*hcSetRealArrXYZ;
		_hc_TYPE_PREFIX T_hcSetRealVecXYZ	*hcSetRealVecXYZ;
		_hc_TYPE_PREFIX T_hcSetStr			*hcSetStr;
		_hc_TYPE_PREFIX T_hcSetStrVecElm		*hcSetStrVecElm;
		_hc_TYPE_PREFIX T_hcSetStrArrElm		*hcSetStrArrElm;
		_hc_TYPE_PREFIX T_hcSetBlock			*hcSetBlock;


		_hc_TYPE_PREFIX T_hcNotifyStart		*hcNotifyStart;
		_hc_TYPE_PREFIX T_hcNotifyStop		*hcNotifyStop;
		_hc_TYPE_PREFIX T_hcNotifySetup		*hcNotifySetup;
		_hc_TYPE_PREFIX T_hcGetNotifyData	*hcGetNotifyData;
		_hc_TYPE_PREFIX T_hcNotifyDataAvail *hcNotifyDataAvail;
		
		_hc_TYPE_PREFIX T_hcAlloc 			*hcAlloc;
		_hc_TYPE_PREFIX T_hcFree 			*hcFree;
		_hc_TYPE_PREFIX T_hcShowMessage 	*hcShowMessage;
		_hc_TYPE_PREFIX T_hcSetTimeouts     *hcSetTimeouts;

		_hc_TYPE_PREFIX T_hcLastError 			*hcLastError;
		_hc_TYPE_PREFIX T_hcGetErrorAction 		*hcGetErrorAction;
		_hc_TYPE_PREFIX T_hcSetErrorAction 		*hcSetErrorAction;

#endif

/***********************************************************************
 * Fortran interface to CDK. So far Fortran stuff do not require	   *
 * typedefs 														   *
 ***********************************************************************/

/**************************************************************************
 * the following CDK routines do not require any change to be called from *
 * FORTRAN.																  *
 **************************************************************************/

	BOOL __stdcall hfConnect(LPSTR lszCmd,unsigned int length_arg);

	BOOL __stdcall hfExecTxt(LPSTR script_cmd,unsigned int length_arg);
	BOOL __stdcall hfQueryTxt(LPSTR var_name,unsigned int length_var,LPSTR res, unsigned int length_res);
	BOOL __stdcall hfQueryBin(HSV hsv,DWORD lwIndx1, DWORD lwIndx2,DWORD* result, DWORD* cbL);
	
	BOOL __stdcall hfNotifyStart(LPSTR var_name,unsigned int length_arg);
	BOOL __stdcall hfNotifyStop(LPSTR var_name,unsigned int length_arg);

	BOOL __stdcall hfStrC2F(LPSTR cS,LPSTR fS,unsigned int length_res);	
	BOOL __stdcall hfStrF2C(LPSTR fS,unsigned int length_arg,LPSTR cS);
	void __stdcall hfShowMessage(LPSTR fs,unsigned int length_arg);
	

#define HC_HEADER_INCLUDED

#endif
