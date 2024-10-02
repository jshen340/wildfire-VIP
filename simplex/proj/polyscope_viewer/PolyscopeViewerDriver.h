//////////////////////////////////////////////////////////////////////////
// Polyscope viewer driver
// Copyright (c) (2018-), Bo Zhu, Zhecheng Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
////This file contains customized viewers for different applications
//////////////////////////////////////////////////////////////////////////

#ifndef __PolyscopeViewerDriver_h__
#define __PolyscopeViewerDriver_h__

#include "PolyscopeViewer.h"
#include "PolyscopeGrid.h"

//////////////////////////////////////////////////////////////////////////
////fluid viewers

class PolyscopeViewerEulerianFluid : public PolyscopeViewer
{typedef PolyscopeViewer Base;
public:
	virtual void Register_Common_Data();
	virtual void Register_Data();
	virtual void Update_Data();
	virtual void Register_UI();
};

#endif
