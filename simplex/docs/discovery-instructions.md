## Request an Account
Go to [https://rc.dartmouth.edu/](https://rc.dartmouth.edu/) and click "Request an Account".
They will send you an e-mail in 1~2 day(s).

## Setup Host Alias
Open `~/.ssh/config`

## Try ssh login
Suppose your NetID is F000XYZ.

        $ ssh discovery

## Install ssh-copy-id
## Copy public key

        $ ssh-copy-id discovery
## Login without password

        $ ssh discovery
## Generate SSH key for GitLab

        $ ssh-keygen -t ed25519 -C "discovery"

## Copy SSH key to GitLab

        $ cat ~/.ssh/id_ed25519.pub

## Clone simplex and complex

        cd ~
        git clone git@gitlab.com:boolzhu/simplex.git
        git clone git@gitlab.com:boolzhu/complex.git
