# Running Cervus using wine

_In theory_ this is quite straight forward. 

Linux users: Install wine using info from https://www.winehq.org/ (confirmed working with ≥ v6.1)

MacOS users: Sorry this gets significantly more hit and miss, especially if you are using an arm64 (Apple silicon) based Mac. For the best chance of success, I recommend the following:

Install the [crossover-wine](https://github.com/Gcenx/winecx) fork of wine with the [winetricks](https://github.com/Gcenx/homebrew-wine) package from github user Gcenx. For a homebrew based guide, see: https://github.com/Gcenx/homebrew-wine. This is probably the easiest way to get things working. If you use an Intel based Mac that can run 32bit programs, then installing wine as normal will probably work fine too. Failing that, the paid version of [Crossover](https://www.codeweavers.com/crossover) is a very convenient method of running wine on MacOS (but obviously has a cost associated with it).


