# Tensorflow Notes

## CPU only

pip install --ignore-installed --upgrade \
https://storage.googleapis.com/tensorflow/linux/cpu/tensorflow-1.6.0-cp36-cp36m-linux_x86_64.whl

## GPU Preparation

### Verify You Have a CUDA-Capable GPU

```bash
$ lspci | grep -i nvidia
04:00.0 VGA compatible controller: NVIDIA Corporation GM204 [GeForce GTX 970] (rev a1)
04:00.1 Audio device: NVIDIA Corporation GM204 High Definition Audio Controller (rev a1)
```

### Verify You have a supported version of Linux

```bash
uname -m && cat /etc/*release
x86_64
DISTRIB_ID=Ubuntu
DISTRIB_RELEASE=17.10
DISTRIB_CODENAME=artful
DISTRIB_DESCRIPTION="Ubuntu 17.10"
NAME="Ubuntu"
VERSION="17.10 (Artful Aardvark)"
ID=ubuntu
ID_LIKE=debian
PRETTY_NAME="Ubuntu 17.10"
VERSION_ID="17.10"
HOME_URL="https://www.ubuntu.com/"
SUPPORT_URL="https://help.ubuntu.com/"
BUG_REPORT_URL="https://bugs.launchpad.net/ubuntu/"
PRIVACY_POLICY_URL="https://www.ubuntu.com/legal/terms-and-policies/privacy-policy"
VERSION_CODENAME=artful
UBUNTU_CODENAME=artful
```


### Verify the System Has gcc Installed

```bash
uname -r
4.13.0-36-generic
```


### Verify the System has the Correct Kernel Headers and Development Packages Installed


```bash
sudo apt-get install linux-headers-$(uname -r)
```


### Install CUDA

```
sudo dpkg -i cuda-repo-ubuntu1704-9-1-local_9.1.85-1_amd64.deb
sudo apt-key add /var/cuda-repo-<version>/7fa2af80.pub
sudo apt-get update
sudo apt-get install cuda
```
