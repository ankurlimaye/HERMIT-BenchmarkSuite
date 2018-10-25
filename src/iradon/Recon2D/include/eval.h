/* Eval.h */ 
extern void KorrelateTest(Image*, Image*, int, int);
extern float LargestValue(Image *MyImage);
extern float SmallestValue(Image *MyImage);
extern float MeanValue(Image *MyImage);
extern float Deviation(Image *MyImage);
extern void PrintStats(int DebugLevel, Image *MyImage);
extern float L1Norm(Image *MyImage1, Image *MyImage2);
extern float L2Norm(Image *MyImage1, Image *MyImage2);
extern Image * DiffImage(Image *MyImage1,Image *MyImage2);

