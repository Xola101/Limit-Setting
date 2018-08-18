//
//   @file    AtlasStyle.h         
//   
//            ATLAS Style, based on a style file from BaBar
//
//
//   @author M.Sutton
// 
//   Copyright (C) 2010 Atlas Collaboration
//
//   $Id: AtlasStyle.h 587003 2014-03-10 17:54:09Z adye $

#ifndef  __ATLASSTYLE_H
#define __ATLASSTYLE_H

class TStyle;

void SetAtlasStyle (bool force=true);

TStyle* AtlasStyle(const char* baseStyle="Modern");

#endif // __ATLASSTYLE_H
