Clone a local copy of the repository as described on the project page and change directory to your new local version of HCV-sequence-capture-info-flow. Check out a new branch with name as e.g., your initials and your branch count.
```
git checkout -b dw_1
```
Edit `HCV_seqcap_info_flow.xml` using the [desktop version](https://github.com/jgraph/drawio-desktop/releases) of draw.io. Save changes to the same file. Press ctrl-a to select all (so you can crop the export around objects without empty spaces), then export to a png file overwriting `HCV_seqcap_info_flow.png`. Optionally do the same for html (not sure how useful this is yet).
```
git add HCV_seqcap_info_flow.xml
git commit -m "update diagram"
git add HCV_seqcap_info_flow.png
git commit -m "update png export of diagram"
```
Push the branch to PHE GitLab.
```
git push --set-upstream origin dw_1
```
Follow the link for the merge request. Follow the instructions and click the remove source branch to do the merge. Checkout the main branch and pull the latest version.
```
git checkout master
git pull
```
Delete your local branch if you want to.
```
git branch rm dw_1
```
