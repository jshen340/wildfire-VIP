set output=%1
set image=%3
set file=%2

if "%output%"=="" goto :error1
if "%file%"=="" set file=out
rem if "%image%"=="" set image=_image

echo "Create simplex animation!"
echo "output folder" %output%"/"%image%
echo "output file name" %file%

ffmpeg -framerate 25 -i %output%/%image%/%%4d.png -c:v libx264 -vf "crop=trunc(iw/2)*2:trunc(ih/2)*2" -pix_fmt yuv420p %output%/%image%/%file%.mp4

goto :finish_ffmpeg_successfully

:error1
echo [Error]: output folder not specified!

:finish_ffmpeg_successfully