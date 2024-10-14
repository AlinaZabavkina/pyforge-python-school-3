# Application Deployment via GitHub Actions

This repository automates the deployment of my python web application to an AWS EC2 instance using GitHub workflow and SSH

## Deployment Workflow
The workflow is triggered by a push to the main or feat/AWS_workflow branch(the second one just for testing purposes).
It performs the following steps:

1. AWS Configuration: The workflow uses AWS IAM user credentials to access AWS resources. The AWS credentials are stored in Github secrets, obviously can't be stored in repo by security reasons.
2. EC2 Deployment: Code is uploaded to an EC2 instance via git clone/git pull. The workflow connects to the instance via SSH using https://github.com/appleboy/ssh-action, installs required software and runs the application.

## How to Deploy

1. Push changes to the main or feat/AWS_workflow branch.
2. Check the logs of 'deploy' action in Github actions. The workflow will automatically deploy the application to the EC2 instance.