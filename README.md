# AMPA-project



## Getting started

To make it easy for you to get started with GitLab, here's a list of recommended next steps.

Already a pro? Just edit this README.md and make it your own. Want to make it easy? [Use the template at the bottom](#editing-this-readme)!

## Add your files

- [ ] [Create](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#create-a-file) or [upload](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#upload-a-file) files
- [ ] [Add files using the command line](https://docs.gitlab.com/ee/gitlab-basics/add-file.html#add-a-file-using-the-command-line) or push an existing Git repository with the following command:

```
cd existing_repo
git remote add origin https://gitlab.rlp.net/surwagle/ampa-project.git
git branch -M main
git push -uf origin main
```

## Integrate with your tools

- [ ] [Set up project integrations](https://gitlab.rlp.net/surwagle/ampa-project/-/settings/integrations)

## Collaborate with your team

- [ ] [Invite team members and collaborators](https://docs.gitlab.com/ee/user/project/members/)
- [ ] [Create a new merge request](https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html)
- [ ] [Automatically close issues from merge requests](https://docs.gitlab.com/ee/user/project/issues/managing_issues.html#closing-issues-automatically)
- [ ] [Enable merge request approvals](https://docs.gitlab.com/ee/user/project/merge_requests/approvals/)
- [ ] [Automatically merge when pipeline succeeds](https://docs.gitlab.com/ee/user/project/merge_requests/merge_when_pipeline_succeeds.html)

## Test and Deploy

Use the built-in continuous integration in GitLab.

- [ ] [Get started with GitLab CI/CD](https://docs.gitlab.com/ee/ci/quick_start/index.html)
- [ ] [Analyze your code for known vulnerabilities with Static Application Security Testing(SAST)](https://docs.gitlab.com/ee/user/application_security/sast/)
- [ ] [Deploy to Kubernetes, Amazon EC2, or Amazon ECS using Auto Deploy](https://docs.gitlab.com/ee/topics/autodevops/requirements.html)
- [ ] [Use pull-based deployments for improved Kubernetes management](https://docs.gitlab.com/ee/user/clusters/agent/)
- [ ] [Set up protected environments](https://docs.gitlab.com/ee/ci/environments/protected_environments.html)

***

# Editing this README

When you're ready to make this README your own, just edit this file and use the handy template below (or feel free to structure it however you want - this is just a starting point!). Thank you to [makeareadme.com](https://www.makeareadme.com/) for this template.

## Suggestions for a good README
Every project is different, so consider which of these sections apply to yours. The sections used in the template are suggestions for most open source projects. Also keep in mind that while a README can be too long and detailed, too long is better than too short. If you think your README is too long, consider utilizing another form of documentation rather than cutting out information.

## Name
Choose a self-explaining name for your project.


# Project

This project contains the scripts used to analyse mRNA localization and generate figures for the paper Wagle et al. 

## Pre-requisite

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install following packages.

```bash
asteval==0.9.28
contourpy==1.0.6
cycler==0.11.0
fonttools==4.38.0
future==0.18.2
kiwisolver==1.4.4
lmfit==1.0.3
matplotlib==3.6.2
numpy==1.23.4
packaging==21.3
pandas==1.5.1
patsy==0.5.3
Pillow==9.3.0
pyparsing==3.0.9
python-dateutil==2.8.2
pytz==2022.6
scikit-posthocs==0.7.0
scipy==1.9.3
seaborn==0.12.1
six==1.16.0
statsmodels==0.13.5
uncertainties==3.1.7
```
To install the above-listed packages, run: 
```bash
pip install -r requirment.txt 
```
Or you can choose to install them separately
## Usage

```python
python mRNA_Analysis_v3.py -m list of mRNAs -w list of widths
```

## Contributing

This code was developed by [Surbhit Wagle](https://sites.google.com/view/surbhitwagle/home)

## License

[MIT](https://choosealicense.com/licenses/mit/)
## Support
This work is supported by the [CRC 1080](https://www.crc1080.com) funded by DFG .


## Authors and acknowledgment
We acknowledge our collaborators and members of [Tchumatchenko Lab](http://tchumatchenko.de)

## License
For open source projects, say how it is licensed.

## Project status
We are curretly still working on updating the scripts for small changes in the plots and any new analysis that might come up.
