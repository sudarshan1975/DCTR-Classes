using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DCTRClasses
{
    public class StringFunctions
    {
        public static List<string>? Split(string inpStr, List<string> splitStrings, out List<string> splitters)
        {
            splitters = [];

            if (inpStr == null) return null;

            List<string> outList = [];

            string curStr = "";
            int i = 0;
            while (i < inpStr.Length)
            {
                curStr += inpStr[i];

                foreach (string s in splitStrings)
                {
                    if (curStr.EndsWith(s))
                    {
                        splitters.Add(s);

                        outList.Add(curStr[..^s.Length]);

                        curStr = "";

                        break;
                    }
                }

                i++;
            }

            outList.Add(curStr);

            return outList;
        }

        public static List<string> SplitString(string inpStr, List<int> countList)
        {
            List<string> outList = [];

            string curStr = "";
            int curCountIx = 0, targetCount = int.MaxValue;
            if (curCountIx < countList.Count) targetCount = countList[curCountIx];

            for (int i = 0; i < inpStr.Length; i++)
            {
                char c = inpStr[i];

                if (targetCount > i) curStr += c;
                else
                {
                    outList.Add(curStr);
                    curStr = "" + c;
                    curCountIx++;
                    if (curCountIx < countList.Count) targetCount += countList[curCountIx];
                    else targetCount = int.MaxValue;
                }
            }

            outList.Add(curStr);

            return outList;
        }
    }
}
