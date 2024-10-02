import os
import sys
import platform

if platform.system().lower() == 'windows':
    sysid='win'
else:
    sysid='unix'

origin_root=os.path.join('..')

if __name__=='__main__':
    if len(sys.argv)==2:
        plex_name,proj_name="simplex",sys.argv[1]
    else:
        assert (len(sys.argv)==3),'Error: must pass 1 or 2 parameter(s). Example: python make_project.py fluid_euler ; python make_project.py simplex fluid_euler'
        plex_name,proj_name=sys.argv[1],sys.argv[2]
        assert (plex_name=='simplex' or plex_name=='complex' or plex_name=='gpu_multigrid'), 'Error: the first parameter must be \"simplex\" or \"complex\"'
    proj_root=os.path.join(origin_root,plex_name,'proj',proj_name)
    build_root=os.path.join(origin_root,plex_name,'build',proj_name)
    if sysid=="win":
        cmd="cmake -S {} -B {} -G \"Visual Studio 17 2022\"".format(proj_root, build_root)
    else:
        cmd="cmake -S {} -B {} -DCMAKE_BUILD_TYPE=Release".format(proj_root, build_root)
    print(cmd)
    os.system(cmd)