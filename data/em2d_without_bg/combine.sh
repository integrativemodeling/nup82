rm -rf combined.png

convert -bordercolor '#ffffff' -border 1x1 no_bg.0.png no_bg.1.png no_bg.2.png no_bg.3.png no_bg.4.png no_bg.5.png no_bg.6.png no_bg.7.png no_bg.8.png no_bg.9.png +append combined1.png

convert -bordercolor '#ffffff' -border 1x1 no_bg.10.png no_bg.12.png no_bg.13.png no_bg.14.png no_bg.15.png no_bg.16.png no_bg.17.png no_bg.18.png no_bg.20.png no_bg.21.png +append combined2.png

convert -bordercolor '#ffffff' -border 1x1 combined1.png combined2.png no_bg.22.png -append combined.png

rm -rf combined1.png
rm -rf combined2.png

