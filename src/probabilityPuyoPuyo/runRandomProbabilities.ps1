# Define constants for steps and L for each dimension
$steps = 1024 * 2
$L2D = 128 * 2
$N_sims = 10

# Define the range of N_species values
$N_species_values = @(6, 7, 8, 9)

# Define the directory containing the executables
$scriptDir = "C:\GitHub\Sequential-Spatial-Killing\src\probabilityPuyoPuyo"

# Get the total number of logical processors and define the maximum number of concurrent jobs
$totalProcessors = (Get-WmiObject -Class Win32_ComputerSystem).NumberOfLogicalProcessors
$maxJobs = [math]::Max(1, $totalProcessors - 4)  # Leave 4 processors free

# Helper function to launch jobs for a given dimension
function Start-DimensionJobs {
    param (
        [string]$exeName,
        [int]$L,
        [int]$steps,
        [array]$N_species_values
    )
    foreach ($N_species in $N_species_values) {
        # Wait until there is a free slot for another job
        while ((Get-Job | Where-Object { $_.State -eq 'Running' }).Count -ge $maxJobs) {
            Start-Sleep -Seconds 1
        }

        # Start a new job for each N_species value
        Start-Job -ScriptBlock {
            param($scriptDir, $exeName, $L, $N_species, $steps)
            & "$scriptDir\$exeName" $L $N_species $steps
            if ($LASTEXITCODE -ne 0) {
                Write-Error "Execution failed for $exeName N_species=$N_species"
            }
        } -ArgumentList $scriptDir, $exeName, $L, $N_species, $steps
    }
}

# Run 2D
for ($sim = 1; $sim -le $N_sims; $sim++) {
    Start-DimensionJobs -exeName "randomProbabilities2D.exe" -L $L2D -steps $steps -N_species_values $N_species_values
}

# Wait for all jobs to complete
Get-Job | Wait-Job

# Clean up the job list
Get-Job | Remove-Job
